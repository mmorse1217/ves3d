template<typename VecContainer>
void ShearFlow(const VecContainer &pos, typename VecContainer::value_type
    shear_rate, VecContainer &vel_inf)
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
            ,stride * sizeof(typename VecContainer::value_type), 
            MemcpyDeviceToDevice);
    }
    axpy(shear_rate, vel_inf, vel_inf); 
}

////////////////////////////////////////////////////////////////////////////////

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
InterfacialVelocity(SurfContainer &S_in, Interaction &Inter, 
    OperatorsMats<value_type, device_type> &mats, 
    const Parameters<value_type> &params, BackgroundFlow &bgFlow) :
    S_(S_in),
    interaction_(Inter),
    bg_flow_(bgFlow),
    Intfcl_force_(params),
    params_(params),
    dt_(params_.ts),
    checked_out_work_sca_(0),
    checked_out_work_vec_(0)
{
    w_sph_.replicate(S_.getPosition());
    velocity.replicate(S_.getPosition());
    tension_.replicate(S_.getPosition());
    
    //Setting initial tension to zero
    tension_.getDevice().Memset(tension_.begin(), 0, 
        tension_.size() * sizeof(value_type));

    int p = S_.getPosition().getShOrder();
    int np = S_.getPosition().getStride();

    //Rot mats
    rot_mat_.resize(p + 1, 1, make_pair(2 * p,np));//to match the rot_chunck
    all_rot_mats_.resize(p + 1, 1, make_pair(np,np));
    all_rot_mats_.getDevice().Memcpy(all_rot_mats_.begin(), mats.all_rot_mats_,
        all_rot_mats_.size() * sizeof(value_type), MemcpyDeviceToDevice);
    
    //W_spherical
    w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
        np * sizeof(value_type), MemcpyDeviceToDevice);
    
    for(int ii=1;ii<S_.getPosition().getNumSubs();++ii)
        w_sph_.getDevice().Memcpy(w_sph_.begin() + ii*np, 
            w_sph_.begin(), np * sizeof(value_type), 
            MemcpyDeviceToDevice);
    
    //Singular quadrature weights
    sing_quad_weights_.resize(1,p);
    sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(), 
        mats.sing_quad_weights_, sing_quad_weights_.size() * sizeof(value_type), 
        MemcpyDeviceToDevice);

    //quadrature weights
    quad_weights_.resize(1,p);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(), mats.quad_weights_,
        quad_weights_.size() * sizeof(value_type), MemcpyDeviceToDevice);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
~InterfacialVelocity()
{
    assert(!checked_out_work_sca_);
    assert(!checked_out_work_vec_);

    purgeTheWorkSpace();
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
updatePositionExplicit(const value_type &dt)
{
    this->dt_ = dt;
    this->updateInteraction();
    
    //Bending
    Vec_t* u1 = checkoutVec();
    Vec_t* u2 = checkoutVec();
    
    Intfcl_force_.bendingForce(S_, *u1);
    stokes(*u1, *u2);   
    axpy(static_cast<value_type>(1), *u2, velocity, velocity);
   
    //Tension
    getTension(velocity, tension_);
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    stokes(*u1, *u2);
    axpy(static_cast<value_type>(1), *u2, velocity, velocity);

    axpy(dt_, velocity, S_.getPosition(), S_.getPositionModifiable());
    
    recycle(u1);
    recycle(u2);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
updatePositionImplicit(const value_type &dt)
{
    this->dt_ = dt;
    this->updateInteraction();

    Vec_t* u1 = checkoutVec();
    Vec_t* u2 = checkoutVec();
    Vec_t* u3 = checkoutVec();

    //Explicit bending for tension
    Intfcl_force_.bendingForce(S_, *u1);
    stokes(*u1, *u2);   
    axpy(static_cast<value_type>(1), *u2, velocity, *u1);
  
    //Tension
    getTension(*u1, tension_);
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    stokes(*u1, *u2);
    axpy(static_cast<value_type>(1), *u2, velocity, *u1);
    
    axpy(dt_, *u1, S_.getPosition(), *u1);
    u2->getDevice().Memcpy(u2->begin(), S_.getPosition().begin(), 
        S_.getPosition().size() * sizeof(value_type), MemcpyDeviceToDevice);
    
    //Update position
    int max_iter(params_.outer_solver_maxit);
    value_type tol(params_.outer_solver_tol);
    enum BiCGSReturn solver_ret;

    int mIter;
    value_type tt;
    int ii, imax(2);

    COUT("  Position solve\n ------------------------------------\n");
    for ( ii=0; ii<imax; ++ii )
    {
        mIter = max_iter;
        tt = tol;
        
        solver_ret = linear_solver_vec_(*this, *u2, *u1, mIter, tt);
        
        if ( solver_ret == BiCGSSuccess )
            break;
        
        if ( (solver_ret  != BiCGSSuccess && ii==imax-1) || 
            solver_ret == RelresIsNan )
        { 
            CERR(" The position solver did not converge with the error \""
                <<solver_ret<<"\"",endl<<endl,exit(1));
            
        }
    }   

    COUTDEBUG(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<< ii * max_iter + mIter
        <<"\n                 Relres = "<<tt<<endl);

    COUTDEBUG("            True relres = "<<
        ((*this)(*u2, *u3),
            axpy(static_cast<value_type>(-1), *u3, *u1, *u3),
            tt = sqrt(AlgebraicDot(*u3, *u3))/sqrt(AlgebraicDot(*u1,*u1))
         )<<endl);

    COUT(" ------------------------------------"<<endl);


    u2->getDevice().Memcpy(S_.getPositionModifiable().begin(), u2->begin(), 
        S_.getPosition().size() * sizeof(value_type), MemcpyDeviceToDevice);

    recycle(u1);
    recycle(u2);
    recycle(u3);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
updateInteraction() const
{
    //Interfacial forces
    Vec_t* u1 = checkoutVec();
    Vec_t* u2 = checkoutVec();
    Vec_t* u3 = checkoutVec();

    Intfcl_force_.bendingForce(S_, *u1);
    Intfcl_force_.tensileForce(S_, tension_, *u3);
    axpy(static_cast<value_type>(1), *u1, *u3, *u3);
    xv(S_.getAreaElement(), *u3, *u3);
    
    //Incorporating the quadrature weights into the density
    for(int ii=0; ii<u3->getNumSubs(); ++ii)
        for(int jj=0; jj<u3->getTheDim(); ++jj)
            u3->getDevice().xy(quad_weights_.begin(), 
                u3->getSubN(ii) + jj * u3->getStride(), u3->getStride(), 
                u3->getSubN(ii) + jj * u3->getStride());                    

    //Self-interaction, to be subtracted
    velocity.getDevice().DirectStokes(S_.getPosition().begin(), u3->begin(), 
        NULL, velocity.getStride(), velocity.getNumSubs(), 
        S_.getPosition().begin(), 0, velocity.getStride(), velocity.begin());

    //Shuffling points and densities
    ShufflePoints(S_.getPosition(), *u1);
    ShufflePoints(*u3, *u2);
    u3->setPointOrder(PointMajor);

    //Far interactions
    typename Interaction::InteractionReturn status;
    status = interaction_(*u1, *u2, *u3);
    
    //Shuffling to the original order
    ShufflePoints(*u3, *u2);
    u1->setPointOrder(AxisMajor);
    u3->setPointOrder(AxisMajor);

    //Subtracting the self-interaction
    if(status)
        axpy(static_cast<value_type>(0), velocity, velocity);
    else
        axpy(static_cast<value_type>(-1), velocity, *u2, velocity);
    
    //Background flow
    bg_flow_(S_.getPosition(), params_.bg_flow_param, *u2);
    axpy(static_cast<value_type>(1), *u2, velocity, velocity);

    recycle(u1);
    recycle(u2);
    recycle(u3);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::getTension(
    const Vec_t &vel_in, Sca_t &tension) const
{
    Sca_t* rhs = checkoutSca();
    Sca_t* wrk = checkoutSca();

    S_.div(vel_in, *rhs);
    
    axpy(static_cast<typename SurfContainer::value_type>(-1), *rhs, *rhs);
    
    int max_iter(params_.inner_solver_maxit);
    value_type tol(params_.inner_solver_tol);
    enum BiCGSReturn solver_ret;
    ///@todo add the restart option to the Parameters
    
    int mIter;
    value_type tt;
    int ii, imax(4);
    COUT("  Tension solve\n ------------------------------------\n");
    for ( ii=0; ii<imax; ++ii )
    {
        mIter = max_iter;
        tt = tol;
        
        solver_ret = linear_solver_(*this, tension, *rhs, mIter, tt);
        if ( solver_ret == BiCGSSuccess )
            break;
        
        if ( (solver_ret  != BiCGSSuccess && ii==imax-1) 
            || solver_ret == RelresIsNan )
            
            CERR(" The tension solver did not converge with the error \""
                <<solver_ret<<"\"",endl<<endl,exit(1));
    }
    COUTDEBUG(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<< ii * max_iter + mIter
        <<"\n                 Relres = "<<tt<<endl);
    COUTDEBUG("            True relres = "<<
        ((*this)(tension, *wrk),
            axpy(static_cast<value_type>(-1), *wrk, *rhs, *wrk),
            tt = sqrt(AlgebraicDot(*wrk, *wrk))/sqrt(AlgebraicDot(*rhs,*rhs))
         )<<endl);

    COUT(" ------------------------------------"<<endl);
    recycle(wrk);
    recycle(rhs);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::stokes(
    const Vec_t &force, Vec_t &velocity) const
{
    int p  = S_.getPosition().getShOrder();
    int np = S_.getPosition().getStride();
    int nv = S_.getPosition().getNumSubs();
    //int rot_chunck = 2 * p * np;
    
    value_type alpha(1.0), beta(0.0);
    int trg_idx(0);
    
    Sca_t* t1 = checkoutSca();;
    Sca_t* t2 = checkoutSca();;

    xyInv(S_.getAreaElement(), w_sph_, *t1);
    int nvX3= force.getTheDim() * nv;

    Vec_t* v1 = checkoutVec();;
    Vec_t* v2 = checkoutVec();;
    
    for(int ii=0;ii <= p; ++ii)
    {
        for(int jj=0;jj < 2 * p; ++jj)
        {
            CircShift(all_rot_mats_.begin() + ii * np * np, jj * np, rot_mat_);
            
            S_.getPosition().getDevice().gemm("N", "N", &np, &nvX3, &np, 
                &alpha, rot_mat_.begin(), &np, S_.getPosition().begin(), 
                &np, &beta, v1->begin(), &np);
            
            S_.getPosition().getDevice().gemm("N", "N", &np, &nvX3, &np, 
                &alpha, rot_mat_.begin(), &np, force.begin(), &np, 
                &beta, v2->begin(), &np);
            
            S_.getPosition().getDevice().gemm("N", "N", &np, &nv, &np,
                &alpha, rot_mat_.begin(), &np, t1->begin(), 
                &np, &beta, t2->begin(), &np);
            
            xy(*t2, w_sph_, *t2);
            xv(*t2, *v2, *v2);
            
            S_.getPosition().getDevice().DirectStokes(v1->begin(), v2->begin(), 
                sing_quad_weights_.begin(), np, nv, S_.getPosition().begin(), 
                ii*2*p + jj, ii*2*p + jj + 1, velocity.begin());
        }
    }
    recycle(t1);
    recycle(t2);
    recycle(v1);
    recycle(v2);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const
{
    Vec_t* fb = checkoutVec();

    Intfcl_force_.linearBendingForce(S_, x_new, *fb);
    stokes(*fb, time_mat_vec);
    axpy(-dt_, time_mat_vec, x_new, time_mat_vec);
    recycle(fb);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::operator()(
    const Sca_t &tension, Sca_t &div_stokes_fs) const
{
    Vec_t* fs = checkoutVec();
    Vec_t* u = checkoutVec();
    
    Intfcl_force_.tensileForce(S_, tension, *fs);
    stokes(*fs, *u);
    S_.div(*u, div_stokes_fs);

    recycle(fs);
    recycle(u);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::reparam()
{
    value_type ts(params_.rep_ts);
    value_type vel;

    int ii(-1);
    Vec_t* u1 = checkoutVec();
    Vec_t* u2 = checkoutVec();
    Sca_t* wrk = checkoutSca();
 
    COUT("  Reparametrization \n ------------------------------------\n");
    while ( ++ii < params_.rep_maxit )
    {
        S_.getSmoothedShapePosition(*u1);
        axpy(static_cast<value_type>(-1), S_.getPosition(), 
            *u1, *u1);
        
        S_.mapToTangentSpace(*u1);
        
        //Advecting tension
        S_.grad(tension_, *u2);
        GeometricDot(*u2, *u1, *wrk);
        axpy(ts, *wrk, tension_, tension_);

        axpy(ts, *u1, S_.getPosition(), S_.getPositionModifiable());
        
        vel = MaxAbs(*u1);
             
        COUTDEBUG("\n              Iteration = "<<ii
            <<"\n                  |vel| = "<<vel<<endl);
        
        if(vel < params_.rep_tol )
            break;

    }
    COUTDEBUG(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<<ii
        <<"\n                  |vel| = "<<vel
        <<"\n ------------------------------------"<<endl);

    recycle(u1);
    recycle(u2);
    recycle(wrk);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
typename SurfContainer::Sca_t* InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
checkoutSca() const
{
    Sca_t* scp;
    
    if(scalar_work_q_.empty())
        scp = new Sca_t;
    else
    {
        scp = scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    scp->replicate(S_.getPosition());
    ++checked_out_work_sca_;
    return(scp);
}
template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
recycle(Sca_t* scp) const
{
    scalar_work_q_.push(scp);
    --checked_out_work_sca_;
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
typename SurfContainer::Vec_t* InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
checkoutVec() const
{
    Vec_t* vcp;
    
    if(vector_work_q_.empty())
        vcp = new Vec_t;
    else
    {
        vcp = vector_work_q_.front();
        vector_work_q_.pop();
    }
    
    vcp->replicate(S_.getPosition());
    ++checked_out_work_vec_;
    
    return(vcp);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
recycle(Vec_t* vcp) const
{
    vector_work_q_.push(vcp);
    --checked_out_work_vec_;
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
void InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
purgeTheWorkSpace() const
{
    while ( !scalar_work_q_.empty() )
    {
         delete scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    while ( !vector_work_q_.empty() )
    {
        delete vector_work_q_.front();
        vector_work_q_.pop();
    }
}

#ifndef NDEBUG

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkExplicit(Vec_t &Fb, Vec_t &SFb, Vec_t &vel, Sca_t &tension, 
    Vec_t &xnew, value_type tol)
{
    bool res = 
        benchmarkBendingForce(S_.getPosition(), Fb, tol) &&
        benchmarkStokes(Fb, SFb, tol)                    &&
        benchmarkBgFlow(SFb, vel, tol)                   &&
        benchmarkTension(vel, tension, tol)              &&
        benchmarkNewPostitionExplicit(xnew, tol);
    
    COUTDEBUG("\n\n -------------------------------"
        <<" Explicit stepper benchmark with " 
        << S_.getPosition().getNumSubs() << " surface(s) "
        <<((res) ? "*Passed*" : "*Failed*" )
        <<" -------------------------------"<<endl);
    return ( res );
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkImplicit(Sca_t &tension, Vec_t &matvec, Vec_t &xnew, value_type tol)
{
    bool res = 
        benchmarkTensionImplicit(tension, tol)                 &&
        benchmarkMatVecImplicit(S_.getPosition(), matvec, tol) &&
        benchmarkNewPostitionImplicit(xnew, tol);

        
    COUTDEBUG("\n\n -------------------------------"
        <<" Implicit stepper benchmark with " 
        << S_.getPosition().getNumSubs() << " surface(s) "
        <<((res) ? "*Passed*" : "*Failed*" )
        <<" -------------------------------"<<endl);
    return ( res );
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkBendingForce(const Vec_t &x, Vec_t &Fb, value_type tol) const
{
    COUTDEBUG("\n  Bending force benchmark"
        <<"\n ------------------------------------\n");
    
    Vec_t* u1 = checkoutVec();
    Intfcl_force_.linearBendingForce(S_, x, *u1);
    axpy(static_cast<value_type>(-1), Fb, *u1, Fb);
    
    value_type err = MaxAbs(Fb); 
    axpy(static_cast<value_type>(1), *u1, Fb);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);
        
    recycle(u1);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkStokes(const Vec_t &F, Vec_t &SF, value_type tol) const
{
    COUTDEBUG("\n  Singular Stokes benchmark"
        <<"\n ------------------------------------\n");
    
    Vec_t* u1 = checkoutVec();
    stokes(F, *u1);
    axpy(static_cast<value_type>(-1), SF, *u1, SF);
    
    value_type err = MaxAbs(SF); 
    axpy(static_cast<value_type>(1), *u1, SF);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    recycle(u1);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkBgFlow(const Vec_t &SFb, Vec_t &vel, value_type tol) const
{
    COUTDEBUG("\n  Background flow benchmark"
        <<"\n ------------------------------------\n");

    this->updateInteraction();
    axpy(static_cast<value_type>(1), SFb, velocity, velocity);

    axpy(static_cast<value_type>(-1), vel, velocity, vel);
    
    value_type err = MaxAbs(vel); 
    axpy(static_cast<value_type>(1), velocity, vel);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkTension(const Vec_t &vel, Sca_t &tension, value_type tol) const
{
    COUTDEBUG("\n  Tension benchmark"
        <<"\n ------------------------------------\n");
    
    Sca_t* wrk = checkoutSca();
    axpy(static_cast<value_type>(0), *wrk, *wrk);
    getTension(vel, *wrk);

    axpy(static_cast<value_type>(-1), tension, *wrk, tension);
    
    value_type err = MaxAbs(tension); 
    axpy(static_cast<value_type>(1), *wrk, tension);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    recycle(wrk);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkNewPostitionExplicit(Vec_t &xnew, value_type tol)
{
    COUTDEBUG("\n  Explicit update benchmark"
        <<"\n ------------------------------------\n");

    updatePositionExplicit(dt_);

    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    value_type err = MaxAbs(xnew); 
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);

    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkTensionImplicit(Sca_t &tension, value_type tol)
{
    COUTDEBUG("\n  Tension (implicit) benchmark"
        <<"\n ------------------------------------\n");
    
    this->updateInteraction();

    Vec_t* u1 = checkoutVec();
    Vec_t* u2 = checkoutVec();
    Sca_t* scp = checkoutSca();
    
    //Explicit bending for tension
    Intfcl_force_.bendingForce(S_, *u1);
    stokes(*u1, *u2);   
    axpy(static_cast<value_type>(1), *u2, velocity, *u1);
  
    //Tension
    axpy(static_cast<value_type>(0), *scp, *scp);
    getTension(*u1, *scp);

    axpy(static_cast<value_type>(-1), tension, *scp, tension);
    
    value_type err = MaxAbs(tension); 
    axpy(static_cast<value_type>(1), *scp, tension);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    recycle(u1);
    recycle(u2);
    recycle(scp);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkMatVecImplicit(const Vec_t &x, Vec_t &matvec, value_type tol)
{
    COUTDEBUG("\n  Implicit MatVec benchmark"
        <<"\n ------------------------------------\n");
    
    Vec_t* b = checkoutVec();
    this->operator()(x, *b);

    axpy(static_cast<value_type>(-1), matvec, *b, matvec);
    value_type err = MaxAbs(matvec); 
    axpy(static_cast<value_type>(1), *b, matvec);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    recycle(b);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction, typename BackgroundFlow>
bool InterfacialVelocity<SurfContainer, Interaction, BackgroundFlow>::
benchmarkNewPostitionImplicit(Vec_t &xnew, value_type tol)
{
    COUTDEBUG("\n  Implicit update benchmark"
        <<"\n ------------------------------------\n");

    updatePositionImplicit(dt_);
    
    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    value_type err = MaxAbs(xnew); 
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);
    
    return (err < tol);
}

#endif //NDEBUG
