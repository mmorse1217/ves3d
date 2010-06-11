template<typename SurfContainer>
InterfacialVelocity<SurfContainer>::InterfacialVelocity(SurfContainer *S_in) :
    S_(S_in)
{
    w_sph_.replicate(S_->getPosition());
    u1_.replicate(S_->getPosition());
    u2_.replicate(S_->getPosition());
    tension_.replicate(S_->getPosition());

    int p = S_->getPosition().getShOrder();
    int np = S_->getPosition().getStride();

    sing_quad_weights_.resize(1,p);
    rot_mat_.resize(p + 1, 1, make_pair(2 * p,np));//to match the rot_chunck
    all_rot_mats_.resize(p + 1, 1, make_pair(np,np));


    DataIO<float,CPU> IO(S_->getPosition().getDevice(),"",0);
    char fname[300];
    sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p);
    IO.ReadData(fname, all_rot_mats_.size(), all_rot_mats_.begin());
    
    sprintf(fname,"precomputed/w_sph_%u_single.txt",p);
    IO.ReadData(fname, np, w_sph_.begin());
  
    for(int ii=1;ii<S_->getPosition().getNumSubs();++ii)
        w_sph_.getDevice().Memcpy(w_sph_.begin() + ii*np, 
            w_sph_.begin(), np * sizeof(typename SurfContainer::value_type), 
            MemcpyDeviceToDevice);

    sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p);
    IO.ReadData(fname, np, sing_quad_weights_.begin());
}

template<typename SurfContainer>
void InterfacialVelocity<SurfContainer>::operator()(const value_type &t, 
    Vec &velocity) const
{
    ///@todo Add Interaction
    BgFlow(S_->getPosition(), velocity);
    Intfcl_force_.BendingForce(*S_, u1_);
    Stokes(u1_, u2_);
    axpy(1.0, velocity, u2_, velocity);
  
    //Get tension
    GetTension(velocity, tension_);
    Intfcl_force_.TensileForce(*S_, tension_, u1_);
    Stokes(u1_, u2_);
    axpy(1.0, velocity, u2_, velocity);
}

template<typename SurfContainer>
void InterfacialVelocity<SurfContainer>::GetTension(const Vec &vel_in,
    Sca &tension) const
{
    Sca rhs(tension.getNumSubs(), tension.getShOrder());
    S_->div(vel_in, rhs);
    
    int max_iter = 1000;
    value_type tol = 1e-6;
    
    typename Sca::iterator it = tension.begin();
    for ( ;it !=tension.end(); ++it)
        *it = 0;

    BiCGStab<Sca, InterfacialVelocity<SurfContainer> > solver;
    if(solver(*this, tension, rhs, max_iter, tol) !=BiCGSSuccess)
    {
        cerr<<"The tension solver did not converge!"<<endl;
        abort();
    }
    //cout<<max_iter<<"\t"<<tol<<endl;
}

template<typename SurfContainer>
void InterfacialVelocity<SurfContainer>::BgFlow(const Vec &pos,
    Vec &vel_inf) const
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;

    axpy(0.0, pos, vel_inf);
    for(int ii=0;ii<n_surfs;ii++)
    {
        idx = 3 * ii * stride;
        pos.getDevice().Memcpy(vel_inf.begin() + idx, 
            pos.begin() + idx + stride + stride
            ,stride * sizeof(typename SurfContainer::value_type), 
            MemcpyDeviceToDevice);
    }
    axpy(0.1, vel_inf, vel_inf); 
}

template<typename SurfContainer>
void InterfacialVelocity<SurfContainer>::Stokes(const Vec &force,
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
    int nvX3= 3 * nv;

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

template<typename SurfContainer>
void InterfacialVelocity<SurfContainer>::operator()(const Sca &tension, 
    Sca &div_stokes_fs) const
{
    Vec Fs(tension.getNumSubs(), tension.getShOrder());
    Vec u(tension.getNumSubs(), tension.getShOrder());
    
    Intfcl_force_.TensileForce(*S_, tension, Fs);
    Stokes(Fs, u);
    S_->div(u, div_stokes_fs);
}

// template<typename SurfContainer>
// void InterfacialVelocity<SurfContainer>::TensionPrecond::operator()(const 
//     Sca &in, Sca &out) const
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
// }


// template <typename ScalarContainer, typename VectorContainer > 
// void Surface<ScalarContainer, VectorContainer>::Reparam()
// {
//     int iter = 0;
//     T vel = 2*params_.rep_max_vel_;

//     //up-sample
//     x_.getDevice().InterpSh(params_.p_, 3*params_.n_surfs_, x_.begin(), work_arr.begin(), shc.begin(), params_.rep_up_freq_, V10.begin());

//     while(iter++ < params_.rep_iter_max_ && vel > params_.rep_max_vel_ * params_.rep_max_vel_)
//     {
//         UpdateNormal();///@todo shc.begin() may have changed
//         x_.getDevice().Filter(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), alpha_q, work_arr.begin(), shc.begin(), V11.begin());

//         axpy((T) -1.0, V10, V11, V11);
//         GeometricDot(V13, V11, S10);
//         axpb((T) -1.0,S10, (T) 0.0, S10);
//         xvpw(S10, V13, V11, V11);//The correction velocity
//         axpy(params_.rep_ts_, V11, V10, V10);

//         GeometricDot(V11,V11,S10);
//         ///@todo The error should be relative
//         vel = S10.Max();
// #ifndef NDEBUG
//         cout<<" Reparam iteration "<<iter<<", max vel: "<<vel<<endl;
// #endif
//     }
    
//     ///inteprsh is buggy
//     x_.getDevice().InterpSh(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), work_arr.begin(), shc.begin(), params_.p_, V11.begin());
//     x_.getDevice().Memcpy(x_.begin(), V11.begin(), x_.getDataLength(), MemcpyDeviceToDevice);
// }
