template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::
InterfacialVelocity(SurfContainer &S_in, const Interaction &Inter, 
    OperatorsMats<Sca_t> &mats, 
    const Parameters<value_type> &params, const BgFlowBase<Vec_t> &bgFlow) :
    usr_ptr_(NULL),
    S_(S_in),
    interaction_(Inter),
    bg_flow_(bgFlow),
    Intfcl_force_(params),
    params_(params),
    dt_(params_.ts),   
    sht_(mats.p_, mats.mats_p_),
    sht_upsample_(mats.p_up_, mats.mats_p_up_),
    move_pole(mats),
    checked_out_work_sca_(0),
    checked_out_work_vec_(0)
{
    velocity_.replicate(S_.getPosition());
    tension_.replicate(S_.getPosition());
    
    //Setting initial tension to zero
    tension_.getDevice().Memset(tension_.begin(), 0, 
        tension_.size() * sizeof(value_type));

    int p = S_.getPosition().getShOrder();
    int np = S_.getPosition().getStride();
    
    //W_spherical
    w_sph_.resize(1, p);
    w_sph_inv_.resize(1, p);
    w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
        np * sizeof(value_type), MemcpyDeviceToDevice);
    xInv(w_sph_,w_sph_inv_);
    
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

template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::
~InterfacialVelocity()
{
    assert(!checked_out_work_sca_);
    assert(!checked_out_work_vec_);

    purgeTheWorkSpace();
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updatePositionExplicit(const value_type &dt)
{
    this->dt_ = dt;
    this->updateInteraction();

    //Bending
    auto_ptr<Vec_t> u1 = checkoutVec();
    auto_ptr<Vec_t> u2 = checkoutVec();
    
    Intfcl_force_.bendingForce(S_, *u1);
    QC(stokes(*u1, *u2));   
    axpy(static_cast<value_type>(1), *u2, velocity_, velocity_);
   
    //Tension
    QC(getTension(velocity_, tension_));
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    QC(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1), *u2, velocity_, velocity_);

    axpy(dt_, velocity_, S_.getPosition(), S_.getPositionModifiable());
    
    recycle(u1);
    recycle(u2);
    
    return Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updatePositionImplicit(const value_type &dt)
{
    this->dt_ = dt;
    this->updateInteraction();
    
    auto_ptr<Vec_t> u1 = checkoutVec();
    auto_ptr<Vec_t> u2 = checkoutVec();
    auto_ptr<Vec_t> u3 = checkoutVec();
    
    //Explicit bending for tension
    Intfcl_force_.bendingForce(S_, *u1);
    QC(stokes(*u1, *u2));   
    axpy(static_cast<value_type>(1), *u2, velocity_, *u1);
  
    //Tension
    QC(getTension(*u1, tension_));
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    QC(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1), *u2, velocity_, *u1);
    
    axpy(dt_, *u1, S_.getPosition(), *u1);
    u2->getDevice().Memcpy(u2->begin(), S_.getPosition().begin(), 
        S_.getPosition().size() * sizeof(value_type), MemcpyDeviceToDevice);
    
    //Update position
    int iter(params_.position_solver_iter);
    int rsrt(params_.position_solver_restart);
    value_type tol(params_.position_solver_tol);
    enum BiCGSReturn solver_ret;
    Error_t ret_val(Success);

    COUT("  Position solve\n ------------------------------------\n");
    solver_ret = linear_solver_vec_(*this, *u2, *u1, rsrt, iter, tol);
    if ( solver_ret  != BiCGSSuccess )
        ret_val = SolverDiverged;    
    COUT(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<< iter
        <<"\n                 Relres = "<<tol<<endl);
    
    COUTDEBUG("            True relres = "<<
        ((*this)(*u2, *u3),
            axpy(static_cast<value_type>(-1), *u3, *u1, *u3),
            tol = sqrt(AlgebraicDot(*u3, *u3))/sqrt(AlgebraicDot(*u1,*u1))
         )<<endl);
    
    COUT(" ------------------------------------"<<endl);
    
    u2->getDevice().Memcpy(S_.getPositionModifiable().begin(), u2->begin(), 
        S_.getPosition().size() * sizeof(value_type), MemcpyDeviceToDevice);
    
    recycle(u1);
    recycle(u2);
    recycle(u3);

    return ret_val;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateInteraction() const
{
    //Interfacial forces
    auto_ptr<Vec_t> u1 = checkoutVec();
    auto_ptr<Vec_t> u2 = checkoutVec();
    auto_ptr<Vec_t> u3 = checkoutVec();
    velocity_.replicate(S_.getPosition());
    auto_ptr<Vec_t> shc = checkoutVec();
    auto_ptr<Vec_t> wrk = checkoutVec();

    Intfcl_force_.bendingForce(S_, *u1);
    Intfcl_force_.tensileForce(S_, tension_, *u3);
    axpy(static_cast<value_type>(1), *u1, *u3, *u3);
    xv(S_.getAreaElement(), *u3, *u3);
    
    //Incorporating the quadrature weights into the density
    ax<Sca_t>(quad_weights_,*u3, *u3);

    //Self-interaction, to be subtracted
    velocity_.getDevice().DirectStokes(S_.getPosition().begin(), u3->begin(), 
        static_cast<const value_type*>(NULL), velocity_.getStride(), 
        velocity_.getNumSubs(), S_.getPosition().begin(), 0, 
        velocity_.getStride(), velocity_.begin());

    ///@todo add the flag to parameters
    if ( params_.upsample_interaction )
    {
        //upsampling
        int usf(sht_upsample_.getShOrder());
        u1->resize(u1->getNumSubs(), usf);
        u2->resize(u2->getNumSubs(), usf);
        shc->resize(shc->getNumSubs(), usf);
        wrk->resize(wrk->getNumSubs(), usf);
    
        Resample(*u3, sht_, sht_upsample_, *shc, *wrk, *u1);
        ShufflePoints(*u1, *u2);

        u3->resize(u3->getNumSubs(), usf);
        Resample(S_.getPosition(), sht_, sht_upsample_, *shc, *wrk, *u3);
        ShufflePoints(*u3, *u1);
        u3->setPointOrder(PointMajor); 
    }
    else
    {
        //Shuffling points and densities
        ShufflePoints(S_.getPosition(), *u1);
        ShufflePoints(*u3, *u2);
        u3->setPointOrder(PointMajor);
        //u2->setPointOrder(PointMajor);
    }
    //Far interactions
    Error_t status = interaction_(*u1, *u2, *u3, usr_ptr_);

    //Shuffling to the original order
    ShufflePoints(*u3, *u2);
    u1->setPointOrder(AxisMajor);
    u3->setPointOrder(AxisMajor);

    if ( params_.upsample_interaction )
    {
        Resample(*u2, sht_upsample_, sht_, *shc, *wrk, *u1);

        int dsf(sht_.getShOrder());
        u1->resize(u1->getNumSubs(), dsf);
        u2->resize(u2->getNumSubs(), dsf);
        u3->resize(u3->getNumSubs(), dsf);
        shc->resize(shc->getNumSubs(), dsf);
        wrk->resize(wrk->getNumSubs(), dsf);
        
        sht_.lowPassFilter(*u1, *wrk, *shc, *u2);
    }

    //Subtracting the self-interaction
    if ( status == Success )
        axpy(static_cast<value_type>(-1), velocity_, *u2, velocity_); 
    else
        axpy(static_cast<value_type>(0), velocity_, velocity_);
    
    //Background flow
    ///@bug the time should be passed to the BgFlow handle.
    bg_flow_(S_.getPosition(), 0, *u2);
    axpy(static_cast<value_type>(1), *u2, velocity_, velocity_);

    recycle(u1);
    recycle(u2);
    recycle(u3);
    recycle(shc);
    recycle(wrk);

    return Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::getTension(
    const Vec_t &vel_in, Sca_t &tension) const
{
    auto_ptr<Sca_t> rhs = checkoutSca();
    auto_ptr<Sca_t> wrk = checkoutSca();

    S_.div(vel_in, *rhs);
    
    axpy(static_cast<value_type>(-1), *rhs, *rhs);
    
    int iter(params_.tension_solver_iter);
    int rsrt(params_.tension_solver_restart);
    value_type tol(params_.tension_solver_tol);
    enum BiCGSReturn solver_ret;
    Error_t ret_val(Success);
    
    COUT("  Tension solve\n ------------------------------------\n");
    solver_ret = linear_solver_(*this, tension, *rhs, rsrt, iter, tol);
    if ( solver_ret  != BiCGSSuccess )
        ret_val = SolverDiverged;
    COUT(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<< iter
        <<"\n                 Relres = "<<tol<<endl);
    COUTDEBUG("            True relres = "<<
        ((*this)(tension, *wrk),
            axpy(static_cast<value_type>(-1), *wrk, *rhs, *wrk),
            tol = sqrt(AlgebraicDot(*wrk, *wrk))/sqrt(AlgebraicDot(*rhs,*rhs))
         )<<endl);
    COUT(" ------------------------------------"<<endl);
    
    recycle(wrk);
    recycle(rhs);
    
    return ret_val;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::stokes(
    const Vec_t &force, Vec_t &velocity) const
{
    PROFILESTART();

    int imax(S_.getPosition().getGridDim().first);
    int jmax(S_.getPosition().getGridDim().second);
    int np = S_.getPosition().getStride();
    int nv = S_.getPosition().getNumSubs();

    auto_ptr<Sca_t> t1 = checkoutSca();
    auto_ptr<Sca_t> t2 = checkoutSca();
    auto_ptr<Vec_t> v1 = checkoutVec();
    auto_ptr<Vec_t> v2 = checkoutVec();

    ax(w_sph_inv_, S_.getAreaElement(), *t1);
    
    int numinputs = 3;
    const Sca_t* inputs[] = {&S_.getPosition(), &force, t1.get()};
    Sca_t* outputs[] = {v1.get(), v2.get(), t2.get()};
    move_pole.setOperands(inputs, numinputs, params_.singular_stokes);

    for(int ii=0;ii < imax; ++ii)
        for(int jj=0;jj < jmax; ++jj)
        {
            move_pole(ii, jj, outputs);

            ax(w_sph_, *t2, *t2); 
            xv(*t2, *v2, *v2);
            
            PROFILESTART();
            S_.getPosition().getDevice().DirectStokes(v1->begin(), v2->begin(), 
                sing_quad_weights_.begin(), np, nv, S_.getPosition().begin(), 
                ii * jmax + jj, ii * jmax + jj + 1, velocity.begin());
            PROFILEEND("SelfInteraction_",0);
        }
    
    recycle(t1);
    recycle(t2);
    recycle(v1);
    recycle(v2);
    
    PROFILEEND("",0);
    return Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const
{
    auto_ptr<Vec_t> fb = checkoutVec();

    Intfcl_force_.linearBendingForce(S_, x_new, *fb);
    QC(stokes(*fb, time_mat_vec));
    axpy(-dt_, time_mat_vec, x_new, time_mat_vec);
    recycle(fb);

    return Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::operator()(
    const Sca_t &tension, Sca_t &div_stokes_fs) const
{
    auto_ptr<Vec_t> fs = checkoutVec();
    auto_ptr<Vec_t> u = checkoutVec();
    
    Intfcl_force_.tensileForce(S_, tension, *fs);
    QC(stokes(*fs, *u));
    S_.div(*u, div_stokes_fs);

    recycle(fs);
    recycle(u);
    
    return Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::reparam()
{
    value_type ts(params_.rep_ts);
    value_type vel;

    int ii(-1);
    auto_ptr<Vec_t> u1 = checkoutVec();
    auto_ptr<Vec_t> u2 = checkoutVec();
    auto_ptr<Sca_t> wrk = checkoutSca();
 
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
             
        COUT("\n              Iteration = "<<ii
            <<"\n                  |vel| = "<<vel<<endl);
        
        if(vel < params_.rep_tol )
            break;

    }
    COUT(" ------------------------------------"<<endl);
    COUT("       Total iterations = "<<ii
        <<"\n                  |vel| = "<<vel
        <<"\n ------------------------------------"<<endl);

    recycle(u1);
    recycle(u2);
    recycle(wrk);
    
    return Success;
}

template<typename SurfContainer, typename Interaction>
auto_ptr<typename SurfContainer::Sca_t> InterfacialVelocity<SurfContainer, Interaction>::
checkoutSca() const
{
    auto_ptr<Sca_t> scp;
    
    if(scalar_work_q_.empty())
        scp = static_cast<auto_ptr<Sca_t> >(new Sca_t);
    else
    {
        scp = static_cast<auto_ptr<Sca_t> >(scalar_work_q_.front());
        scalar_work_q_.pop();
    }
    
    scp->replicate(S_.getPosition());
    ++checked_out_work_sca_;
    return(scp);
}
template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::
recycle(auto_ptr<Sca_t> scp) const
{
    scalar_work_q_.push(scp.release());
    --checked_out_work_sca_;
}

template<typename SurfContainer, typename Interaction>
auto_ptr<typename SurfContainer::Vec_t> InterfacialVelocity<SurfContainer, Interaction>::
checkoutVec() const
{
    auto_ptr<Vec_t> vcp;
    
    if(vector_work_q_.empty())
        vcp = static_cast<auto_ptr<Vec_t> >(new Vec_t);
    else
    {
        vcp = static_cast<auto_ptr<Vec_t> >(vector_work_q_.front());
        vector_work_q_.pop();
    }
    
    vcp->replicate(S_.getPosition());
    ++checked_out_work_vec_;
    
    return(vcp);
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::
recycle(auto_ptr<Vec_t> vcp) const
{
    vector_work_q_.push(vcp.release());
    --checked_out_work_vec_;
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkExplicit(Vec_t &Fb, Vec_t &SFb, Vec_t &vel, Sca_t &tension, 
    Vec_t &xnew, value_type tol)
{
    bool res = 
        benchmarkBendingForce(S_.getPosition(), Fb, tol) &&
        benchmarkStokes(Fb, SFb, tol)                    &&
        benchmarkBgFlow(SFb, vel, tol)                   &&
        benchmarkTension(vel, tension, tol)              &&
        benchmarkNewPostitionExplicit(xnew, tol);
    
    COUT("\n\n -------------------------------"
        <<" Explicit stepper benchmark with " 
        << S_.getPosition().getNumSubs() << " surface(s) "
        <<((res) ? "*Passed*" : "*Failed*" )
        <<" -------------------------------"<<endl);
    return ( res );
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkImplicit(Sca_t &tension, Vec_t &matvec, Vec_t &xnew, value_type tol)
{
    bool res = 
        benchmarkTensionImplicit(tension, tol)                 &&
        benchmarkMatVecImplicit(S_.getPosition(), matvec, tol) &&
        benchmarkNewPostitionImplicit(xnew, tol);

        
    COUT("\n\n -------------------------------"
        <<" Implicit stepper benchmark with " 
        << S_.getPosition().getNumSubs() << " surface(s) "
        <<((res) ? "*Passed*" : "*Failed*" )
        <<" -------------------------------"<<endl);
    return ( res );
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkBendingForce(const Vec_t &x, Vec_t &Fb, value_type tol) const
{
    COUTDEBUG("\n  Bending force benchmark"
        <<"\n ------------------------------------\n");
    
    auto_ptr<Vec_t> u1 = checkoutVec();
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkStokes(const Vec_t &F, Vec_t &SF, value_type tol) const
{
    COUTDEBUG("\n  Singular Stokes benchmark"
        <<"\n ------------------------------------\n");
    
    auto_ptr<Vec_t> u1 = checkoutVec();
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkBgFlow(const Vec_t &SFb, Vec_t &vel, value_type tol) const
{
    COUTDEBUG("\n  Background flow benchmark"
        <<"\n ------------------------------------\n");

    this->updateInteraction();
    axpy(static_cast<value_type>(1), SFb, velocity_, velocity_);

    axpy(static_cast<value_type>(-1), vel, velocity_, vel);
    
    value_type err = MaxAbs(vel); 
    axpy(static_cast<value_type>(1), velocity_, vel);
    
    COUTDEBUG("  The benchmark " 
        <<((err<tol) ? "*Passed*" : "*Failed*" )
        <<" with\n                  error = "<<PRINTFRMT<<err<<endl);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkTension(const Vec_t &vel, Sca_t &tension, value_type tol) const
{
    COUTDEBUG("\n  Tension benchmark"
        <<"\n ------------------------------------\n");
    
    auto_ptr<Sca_t> wrk = checkoutSca();
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkTensionImplicit(Sca_t &tension, value_type tol)
{
    COUTDEBUG("\n  Tension (implicit) benchmark"
        <<"\n ------------------------------------\n");
    
    this->updateInteraction();

    auto_ptr<Vec_t> u1 = checkoutVec();
    auto_ptr<Vec_t> u2 = checkoutVec();
    auto_ptr<Sca_t> scp = checkoutSca();
    
    //Explicit bending for tension
    Intfcl_force_.bendingForce(S_, *u1);
    stokes(*u1, *u2);   
    axpy(static_cast<value_type>(1), *u2, velocity_, *u1);
  
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkMatVecImplicit(const Vec_t &x, Vec_t &matvec, value_type tol)
{
    COUTDEBUG("\n  Implicit MatVec benchmark"
        <<"\n ------------------------------------\n");
    
    auto_ptr<Vec_t> b = checkoutVec();
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

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
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
