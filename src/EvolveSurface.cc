template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::InterfacialVelocity(SurfContainer *S_in,
                                                        value_type* data) :
    S_(S_in)
{
    w_sph_.replicate(S_->getPosition());
    u1_.replicate(S_->getPosition());
    u2_.replicate(S_->getPosition());
    u3_.replicate(S_->getPosition());
    tension_.replicate(S_->getPosition());

    int p = S_->getPosition().getShOrder();
    int np = S_->getPosition().getStride();

    DataIO<value_type,CPU> IO(S_->getPosition().getDevice(),"",0);

    //Rot mats
    rot_mat_.resize(p + 1, 1, make_pair(2 * p,np));//to match the rot_chunck
    all_rot_mats_.resize(p + 1, 1, make_pair(np,np));
    char fname[300];
    sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p);
    IO.ReadData(fname, all_rot_mats_.size(), all_rot_mats_.begin());
    
    //W_spherical
    sprintf(fname,"precomputed/w_sph_%u_single.txt",p);
    IO.ReadData(fname, np, w_sph_.begin());
  
    for(int ii=1;ii<S_->getPosition().getNumSubs();++ii)
        w_sph_.getDevice().Memcpy(w_sph_.begin() + ii*np, 
            w_sph_.begin(), np * sizeof(value_type), 
            MemcpyDeviceToDevice);

    //Singular quadrature weights
    sing_quad_weights_.resize(1,p);
    sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p);
    IO.ReadData(fname, np, sing_quad_weights_.begin());

    //quadrature weights
    quad_weights_.resize(1,p);
    sprintf(fname,"precomputed/quad_weights_%u_single.txt",p);
    IO.ReadData(fname, np, quad_weights_.begin());
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::operator()(const value_type &t, 
    Vec &velocity) const
{
    //Interaction
    Intfcl_force_.BendingForce(*S_, u1_);
    Intfcl_force_.TensileForce(*S_, tension_, u3_);
    axpy(static_cast<value_type>(1), u1_, u3_, u3_);
    xv(S_->getAreaElement(), u3_, u3_);
    
    //Self-interaction, to be subtracted
    velocity.getDevice().DirectStokes(S_->getPosition().begin(), u3_.begin(), 
        quad_weights_.begin(), velocity.getStride(), velocity.getNumSubs(), 
        S_->getPosition().begin(), 0, velocity.getStride(), velocity.begin());

    //Incorporating the quadrature weights into the density
    for(int ii=0; ii<u3_.getNumSubs(); ++ii)
        for(int jj=0; jj<u3_.getTheDim(); ++jj)
            u3_.getDevice().xy(quad_weights_.begin(), u3_.getSubN(ii) + 
                jj * u3_.getStride(), u3_.getStride(), u3_.begin());                    

    //Shuffling points
    ShufflePoints(S_->getPosition(), u1_);
    ShufflePoints(u3_, u2_);
 
    //Far interactions
    interaction_(u1_, u2_, u3_);
    
    //Shuffling to the original order
    ShufflePoints(u3_, u2_);

    //Subtracting the self-interaction
    axpy(static_cast<value_type>(-1), velocity, u2_, velocity);
    //cout<<velocity<<endl;

    //Background
    BgFlow(S_->getPosition(), u2_);
    axpy(static_cast<value_type>(1), u2_, velocity, velocity);
    
    //Bending
    Intfcl_force_.BendingForce(*S_, u1_);
    Stokes(u1_, u2_);   
    axpy(static_cast<value_type>(1), u2_, velocity, velocity);
  
    //Get tension
    GetTension(velocity, tension_);
    Intfcl_force_.TensileForce(*S_, tension_, u1_);
    Stokes(u1_, u2_);
    axpy(static_cast<value_type>(1), u2_, velocity, velocity);
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::GetTension(const Vec &vel_in,
    Sca &tension) const
{
    Sca rhs;
    rhs.replicate(S_->getPosition());
    S_->div(vel_in, rhs);
  
    axpy(static_cast<typename SurfContainer::value_type>(-1), rhs, rhs);
    
    int max_iter(Parameters<typename SurfContainer::value_type>::
        getInstance().inner_solver_maxit);
    value_type tol(Parameters<typename SurfContainer::value_type>::
        getInstance().inner_solver_tol);

    typename Sca::iterator it = tension.begin();
    for ( ;it !=tension.end(); ++it)
        *it = 0;

    BiCGStab<Sca, InterfacialVelocity<SurfContainer, Interaction> > solver;
    if(solver(*this, tension, rhs, max_iter, tol) !=BiCGSSuccess)
    {
        cerr<<"The tension solver did not converge!"<<endl;
        exit(1);
    }
    cout<<max_iter<<"\t"<<tol<<endl;
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::BgFlow(const Vec &pos,
    Vec &vel_inf) const
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;

    axpy(0.0, pos, vel_inf);
    for(int ii=0;ii<n_surfs;ii++)
    {
        idx = pos.getTheDim() * ii * stride;
        pos.getDevice().Memcpy(vel_inf.begin() + idx, 
            pos.begin() + idx + stride + stride
            ,stride * sizeof(value_type), 
            MemcpyDeviceToDevice);
    }
    axpy(Parameters<value_type>::getInstance().bg_flow_param, vel_inf, vel_inf); 
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::Stokes(const Vec &force,
    Vec &velocity) const
{
    int p  = S_->getPosition().getShOrder();
    int np = S_->getPosition().getStride();
    int nv = S_->getPosition().getNumSubs();
    //int rot_chunck = 2 * p * np;
    
    value_type alpha(1.0), beta(0.0);
    int trg_idx(0);
    
    Sca t1(nv, p);
    Sca t2(nv, p);

    xyInv(S_->getAreaElement(), w_sph_, t1);
    int nvX3= force.getTheDim() * nv;

    Vec v1;
    Vec v2;
    
    v1.replicate(S_->getPosition());
    v2.replicate(S_->getPosition());

    for(int ii=0;ii <= p; ++ii)
    {
        for(int jj=0;jj < 2 * p; ++jj)
        {
            CircShift(all_rot_mats_.begin() + ii * np * np, jj * np, rot_mat_);
            
            S_->getPosition().getDevice().gemm("N", "N", &np, &nvX3, &np, 
                &alpha, rot_mat_.begin(), &np, S_->getPosition().begin(), 
                &np, &beta, v1.begin(), &np);
            
            S_->getPosition().getDevice().gemm("N", "N", &np, &nvX3, &np, 
                &alpha, rot_mat_.begin(), &np, force.begin(), &np, 
                &beta, v2.begin(), &np);
            
            S_->getPosition().getDevice().gemm("N", "N", &np, &nv, &np,
                &alpha, rot_mat_.begin(), &np, t1.begin(), 
                &np, &beta, t2.begin(), &np);
            
            xy(t2, w_sph_, t2);
            xv(t2, v2, v2);
            
            S_->getPosition().getDevice().DirectStokes(v1.begin(), v2.begin(), 
                sing_quad_weights_.begin(), np, nv, S_->getPosition().begin(), 
                ii*2*p + jj, ii*2*p + jj + 1, velocity.begin());
        }
    }
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::operator()(const Sca &tension, 
    Sca &div_stokes_fs) const
{
    Vec Fs(tension.getNumSubs(), tension.getShOrder());
    Vec u(tension.getNumSubs(), tension.getShOrder());
    
    Intfcl_force_.TensileForce(*S_, tension, Fs);
    Stokes(Fs, u);
    S_->div(u, div_stokes_fs);
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::TensionPrecond::operator()(const 
    Sca &in, Sca &out) const
{
    axpy(static_cast<value_type>(1), in, out);
// {
// //     if(isempty(L))
// //         eig = @(n) (2*n.^2+2*n-1).*n.*(n+1)./(4*n.^2-1)./(2*n+3);
// //     d1 = size(vecIn,1);
// //     p = (sqrt(8*d1+1)-3)/4;
// //     lambda = eig((0:p)');
// //   lambda = repmat(lambda,1,2*p+1);
// //   lambda(abs(lambda)<1e-10) = 1;
// //   L = 1./lambda(:);
  
// // %   for n=0:p
// // %     for m=-p:p
// // %       ind(n+1,p+m+1) = abs(m)<=n;
// // %     end
// // %   end
// // %   ind = find(~ind);
// // end

// // vecOut = L.*vecIn;
}
