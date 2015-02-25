template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::
InterfacialVelocity(SurfContainer &S_in, const Interaction &Inter,
    OperatorsMats<Arr_t> &mats,
    const Parameters<value_type> &params, const BgFlowBase<Vec_t> &bgFlow,
    PSolver_t *parallel_solver) :
    S_(S_in),
    interaction_(Inter),
    bg_flow_(bgFlow),
    Intfcl_force_(params),
    params_(params),
    //
    parallel_solver_(parallel_solver),
    psolver_configured_(false),
    parallel_matvec_(NULL),
    parallel_rhs_(NULL),
    parallel_u_(NULL),
    //
    dt_(params_.ts),
    sht_(mats.p_, mats.mats_p_),
    sht_upsample_(mats.p_up_, mats.mats_p_up_),
    move_pole(mats),
    checked_out_work_sca_(0),
    checked_out_work_vec_(0),
    stokes_(mats,params_)
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
        np * sizeof(value_type), device_type::MemcpyDeviceToDevice);
    xInv(w_sph_,w_sph_inv_);

    //Singular quadrature weights
    sing_quad_weights_.resize(1,p);
    sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
        mats.sing_quad_weights_, sing_quad_weights_.size() *
        sizeof(value_type),
        device_type::MemcpyDeviceToDevice);

    //quadrature weights
    quad_weights_.resize(1,p);
    quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
        mats.quad_weights_,
        quad_weights_.size() * sizeof(value_type),
        device_type::MemcpyDeviceToDevice);

    int p_up = sht_upsample_.getShOrder();
    quad_weights_up_.resize(1, p_up);
    quad_weights_up_.getDevice().Memcpy(quad_weights_up_.begin(),
        mats.quad_weights_p_up_,
        quad_weights_up_.size() * sizeof(value_type),
        device_type::MemcpyDeviceToDevice);
}

template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::
~InterfacialVelocity()
{
    COUTDEBUG("Destroying an instance of interfacial velocity");
    assert(!checked_out_work_sca_);
    assert(!checked_out_work_vec_);

    purgeTheWorkSpace();

    COUTDEBUG("Deleting parallel matvec and containers");
    delete parallel_matvec_;
    delete parallel_rhs_;
    delete parallel_u_;
}

// Performs the following computation:
// velocity(n+1) = updateFarField( bending(n), tension(n) );    // Far-field
// velocity(n+1)+= stokes( bending(n) );                        // Add near-field due to bending
// tension(n+1)  = getTension( velocity(n+1) )                  // Linear solve to compute new tension
// velocity(n+1)+= stokes( tension(n+1) )                       // Add near-field due to tension
// position(n+1) = position(n) + dt*velocity(n+1)
//
// Notes: tension solve is block implicit.
template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateJacobiExplicit(const value_type &dt)
{
    this->dt_ = dt;

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();

    // puts u_inf and interaction in velocity_
    this->updateFarField();

    // add S[f_b]
    Intfcl_force_.bendingForce(S_, *u1);
    CHK(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1.0), *u2, velocity_, velocity_);

    // compute tension
    CHK(getTension(velocity_, tension_));

    // add S[f_sigma]
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    CHK(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1.0), *u2, velocity_, velocity_);

    axpy(dt_, velocity_, S_.getPosition(), S_.getPositionModifiable());

    recycle(u1);
    recycle(u2);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateJacobiGaussSeidel(const value_type &dt)
{
    this->dt_ = dt;

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();
    std::auto_ptr<Vec_t> u3 = checkoutVec();

    // put far field in velocity_ and the sum with S[f_b] in u1
    this->updateFarField();
    Intfcl_force_.bendingForce(S_, *u2);
    CHK(stokes(*u2, *u1));
    axpy(static_cast<value_type>(1.0), velocity_, *u1, *u1);

    // tension
    CHK(getTension(*u1, tension_));
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    CHK(stokes(*u1, *u2));

    // position rhs
    axpy(static_cast<value_type>(1.0), velocity_, *u2, *u1);
    axpy(dt_, *u1, S_.getPosition(), *u1);

    // initial guess
    u2->getDevice().Memcpy(u2->begin(), S_.getPosition().begin(),
        S_.getPosition().size() * sizeof(value_type),
        u2->getDevice().MemcpyDeviceToDevice);

    int iter(params_.position_solver_iter);
    int rsrt(params_.position_solver_restart);
    value_type tol(params_.position_solver_tol),relres(params_.position_solver_tol);

    enum BiCGSReturn solver_ret;
    Error_t ret_val(ErrorEvent::Success);

    COUTDEBUG("Solving for position");
    solver_ret = linear_solver_vec_(*this, *u2, *u1, rsrt, iter, relres);
    if ( solver_ret  != BiCGSSuccess )
        ret_val = ErrorEvent::DivergenceError;

    COUTDEBUG("Position solve: Total iter = "<<iter<<", relres = "<<tol);
    COUTDEBUG("Checking true relres");
    ASSERT(((*this)(*u2, *u3),
            axpy(static_cast<value_type>(-1), *u3, *u1, *u3),
            relres = sqrt(AlgebraicDot(*u3, *u3))/sqrt(AlgebraicDot(*u1,*u1)),
            relres<tol
           ),
           "relres ("<<relres<<")<tol("<<tol<<")"
           );

    u2->getDevice().Memcpy(S_.getPositionModifiable().begin(), u2->begin(),
        S_.getPosition().size() * sizeof(value_type),
        u2->getDevice().MemcpyDeviceToDevice);

    recycle(u1);
    recycle(u2);
    recycle(u3);

    return ret_val;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateImplicit(const value_type &dt)
{
    this->dt_ = dt;
    SolverScheme scheme(GloballyImplicit);
    INFO("Updating using a step with "<<scheme);
    CHK(Prepare(scheme));
    CHK(AssembleRhs(parallel_rhs_, dt_, scheme));
    CHK(AssembleInitial(parallel_u_, dt_, scheme));
    CHK(Solve(parallel_rhs_, parallel_u_, dt_, scheme));
    CHK(Update(parallel_u_));

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::Prepare(const SolverScheme &scheme) const
{
    COUTDEBUG("Resizing the containers");

    if (velocity_.size() != S_.getPosition().size() ||
	dl_coeff_.size() != S_.getPosition().getNumSubs() ){

      velocity_.replicate(S_.getPosition());
      tension_.replicate(S_.getPosition());

      INFO("zeroing content of velocity and tension arrays");
      velocity_.getDevice().Memset(velocity_.begin(), 0, sizeof(value_type)*velocity_.size());
      tension_.getDevice().Memset(tension_.begin(), 0, sizeof(value_type)*tension_.size());

      //Permitting the viscosity constrast to be different for each vesicle
      size_t     nves(velocity_.getNumSubs());
      value_type lambda(params_.viscosity_contrast);
      value_type *buffer = new  value_type[nves];

      dl_coeff_.resize(nves);
      for (size_t iV(0); iV<nves; ++iV) buffer[iV] = (1.0 - lambda);
      dl_coeff_.getDevice().Memcpy(dl_coeff_.begin(),
				   buffer,
				   nves * sizeof(value_type),
				   device_type::MemcpyHostToDevice);

      vel_coeff_.resize(nves);
      for (size_t iV(0); iV<nves; ++iV) buffer[iV] = 0.5*(1.0 + lambda);
      vel_coeff_.getDevice().Memcpy(vel_coeff_.begin(),
				    buffer,
				    nves * sizeof(value_type),
				    device_type::MemcpyHostToDevice);
      delete[] buffer;
    }

    ASSERT(  velocity_.size() ==   S_.getPosition().size(), "inccorrect size");
    ASSERT( 3*tension_.size() ==   S_.getPosition().size(), "inccorrect size");
    ASSERT(  dl_coeff_.size() ==   S_.getPosition().getNumSubs(), "inccorrect size");
    ASSERT( vel_coeff_.size() ==   S_.getPosition().getNumSubs(), "inccorrect size");

    stokes_.SetSrcCoord(S_);
    stokes_.SetTrgCoord(S_);

    if (scheme==GloballyImplicit){
	ASSERT(parallel_solver_ != NULL, "need a working parallel solver");
	if (!psolver_configured_) CHK(ConfigureSolver(scheme));
    }

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ConfigureSolver(const SolverScheme &scheme) const
{
    ASSERT(scheme==GloballyImplicit, "Unsupported scheme");
    COUTDEBUG("Configuring the parallel solver");

    typedef typename PSolver_t::matvec_type POp;
    typedef typename PSolver_t::vec_type PVec;
    typedef typename PVec::size_type size_type;

    // Setting up the operator
    size_t sz(velocity_.size()+tension_.size());
    CHK(parallel_solver_->LinOpFactory(&parallel_matvec_));
    CHK(parallel_matvec_->SetSizes(sz,sz));
    CHK(parallel_matvec_->SetName("Vesicle interaction"));
    CHK(parallel_matvec_->SetContext(static_cast<const void*>(this)));
    CHK(parallel_matvec_->SetApply(ImplicitApply));
    CHK(parallel_matvec_->Configure());

    // setting up the rhs
    CHK(parallel_solver_->VecFactory(&parallel_rhs_));
    CHK(parallel_rhs_->SetSizes(sz));
    CHK(parallel_rhs_->SetName("rhs"));
    CHK(parallel_rhs_->Configure());

    CHK(parallel_rhs_->ReplicateTo(&parallel_u_));
    CHK(parallel_u_->SetName("solution"));

    // setting up the solver
    CHK(parallel_solver_->SetOperator(parallel_matvec_));
    CHK(parallel_solver_->Configure());

    psolver_configured_ = true;
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
AssembleRhs(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const
{
    ASSERT(scheme==GloballyImplicit, "Unsupported scheme");

    // rhs=[u_inf+Bx;0]
    COUTDEBUG("Evaluate background flow");
    std::auto_ptr<Vec_t> vRhs = checkoutVec();
    vRhs->replicate(S_.getPosition());
    CHK(BgFlow(*vRhs, dt));

    COUTDEBUG("Computing the far-field interaction due to explicit traction jump");
    std::auto_ptr<Vec_t> f  = checkoutVec();
    std::auto_ptr<Vec_t> Sf = checkoutVec();
    Intfcl_force_.bendingForce(S_, *f);
    stokes_.SetDensitySL(f.get());
    stokes_.SetDensityDL(NULL);
    stokes_(*Sf);
    axpy(static_cast<value_type>(1.0), *Sf, *vRhs, *vRhs);

    COUTDEBUG("Computing rhs for div(u)");
    std::auto_ptr<Sca_t> tRhs = checkoutSca();
    tRhs->getDevice().Memset(tRhs->begin(), 0, tRhs->size() * sizeof(value_type));

    ASSERT( vRhs->getDevice().isNumeric(vRhs->begin(), vRhs->size()), "Non-numeric rhs");
    ASSERT( tRhs->getDevice().isNumeric(tRhs->begin(), tRhs->size()), "Non-numeric rhs");

    COUTDEBUG("Copy data to parallel rhs array");
    size_t xsz(vRhs->size()), tsz(tRhs->size());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;
    CHK(parallel_rhs_->GetArray(i, rsz));
    ASSERT(rsz==xsz+tsz,"Bad sizes");
    vRhs->getDevice().Memcpy(i    , vRhs->begin(), xsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    tRhs->getDevice().Memcpy(i+xsz, tRhs->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    CHK(parallel_rhs_->RestoreArray(i));

    recycle(vRhs);
    recycle(f);
    recycle(Sf);
    recycle(tRhs);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
AssembleInitial(PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const
{
    COUTDEBUG("Using current position/tension as initial guess");
    size_t vsz(velocity_.size()), tsz(tension_.size());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;

    ASSERT( velocity_.getDevice().isNumeric(velocity_.begin(), velocity_.size()), "Non-numeric velocity");
    ASSERT( tension_.getDevice().isNumeric(tension_.begin(), tension_.size()), "Non-numeric tension");

    CHK(parallel_u_->GetArray(i, rsz));
    ASSERT(rsz==vsz+tsz,"Bad sizes");
    COUTDEBUG("Copy data to parallel solution array");
    velocity_.getDevice().Memcpy(   i, velocity_.begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    tension_.getDevice().Memcpy(i+vsz, tension_.begin() , tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    CHK(parallel_u_->RestoreArray(i));

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ImplicitApply(const POp_t *o, const value_type *x, value_type *y)
{
    const InterfacialVelocity *F(NULL);
    o->Context((const void**) &F);
    size_t vsz(F->velocity_.size()), tsz(F->tension_.size());

    std::auto_ptr<Vec_t> vel = F->checkoutVec();
    std::auto_ptr<Sca_t> ten = F->checkoutSca();
    std::auto_ptr<Vec_t> fb  = F->checkoutVec();
    std::auto_ptr<Vec_t> fs  = F->checkoutVec();
    std::auto_ptr<Vec_t> Sf  = F->checkoutVec();
    std::auto_ptr<Sca_t> Dv  = F->checkoutSca();
    vel->replicate(F->velocity_);
    ten->replicate(F->tension_);
    fb->replicate(*vel);
    fs->replicate(*vel);
    Sf->replicate(*vel);
    Dv->replicate(*vel);

    COUTDEBUG("Unpacking the input parallel vector");
    vel->getDevice().Memcpy(vel->begin(), x    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    ten->getDevice().Memcpy(ten->begin(), x+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);

    COUTDEBUG("Computing the div term");
    F->S_.div(*vel, *Dv);

    COUTDEBUG("Computing the interfacial forces");
    F->Intfcl_force_.linearBendingForce(F->S_, *vel, *fb);
    F->Intfcl_force_.tensileForce(F->S_, *ten, *fs);
    axpy(static_cast<value_type>(F->dt_), *fb, *fs, *fb);
    F->stokes_.SetDensitySL(fb.get());

    if( fabs(F->params_.viscosity_contrast-1.0)>1e-12){
      COUTDEBUG("Setting the double layer density");
      av(F->dl_coeff_, *vel, *fs);
      F->stokes_.SetDensityDL(fs.get());
    };

    F->stokes_(*Sf);
    if( fabs(F->params_.viscosity_contrast-1.0)>1e-12)
      av(F->vel_coeff_, *vel, *vel);

    axpy(static_cast<value_type>(-1.0), *Sf, *vel, *vel);

    ASSERT(vel->getDevice().isNumeric(vel->begin(), vel->size()), "Non-numeric velocity");
    ASSERT(Dv->getDevice().isNumeric(Dv->begin(), Dv->size()), "Non-numeric divergence");

    vel->getDevice().Memcpy(y   , vel->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    Dv->getDevice().Memcpy(y+vsz, Dv->begin() , tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);

    F->recycle(vel);
    F->recycle(ten);
    F->recycle(fb);
    F->recycle(fs);
    F->recycle(Sf);
    F->recycle(Dv);

    return ErrorEvent::Success;
}


template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
Solve(const PVec_t *rhs, PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const
{
    INFO("Solving for position and tension using "<<scheme<<" scheme.");
    //parallel_rhs_->View();
    //parallel_u_->View();

    CHK(parallel_solver_->Solve(parallel_rhs_, parallel_u_));
    typename PVec_t::size_type	iter;
    CHK(parallel_solver_->IterationNumber(iter));
    INFO("Parallel solver returned after "<<iter<<" iteration(s).");
    INFO("Parallal solver report:");
    parallel_solver_->ViewReport();
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::Update(PVec_t *u0)
{
    COUTDEBUG("Updating position and tension.");
    size_t  vsz(velocity_.size()), tsz(tension_.size());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;

    CHK(parallel_u_->GetArray(i, rsz));
    ASSERT(rsz==vsz+tsz,"Bad sizes");
    COUTDEBUG("Copy data from parallel solution array.");
    velocity_.getDevice().Memcpy(velocity_.begin(), i    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    tension_.getDevice().Memcpy(tension_.begin()  , i+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    CHK(parallel_u_->RestoreArray(i));
    axpy(dt_, velocity_, S_.getPosition(), S_.getPositionModifiable());

    return ErrorEvent::Success;
}

template<typename	SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
BgFlow(Vec_t &bg, const value_type &dt) const{
    //!@bug the time should be passed to the BgFlow handle.
    bg_flow_(S_.getPosition(), 0, bg);

    return ErrorEvent::Success;
}

// Compute velocity_far = velocity_bg + FMM(bending+tension) - DirectStokes(bending+tension)
template<typename	SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateFarField() const
{
    velocity_.replicate(S_.getPosition());
    CHK(this->BgFlow(velocity_, this->dt_));

    if (this->interaction_.HasInteraction()){
	std::auto_ptr<Vec_t>	fi  = checkoutVec();
	std::auto_ptr<Vec_t>	vel = checkoutVec();
	fi->replicate(velocity_);
	vel->replicate(velocity_);

	Intfcl_force_.bendingForce(S_, *fi);
	Intfcl_force_.tensileForce(S_, tension_, *vel);
	axpy(static_cast<value_type>(1.0), *fi, *vel, *fi);

	EvaluateFarInteraction(S_.getPosition(), *fi, *vel);
	axpy(static_cast<value_type>(1.0), *vel, velocity_, velocity_);

	recycle(fi);
	recycle(vel);
    }

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
CallInteraction(const Vec_t &src, const Vec_t &den, Vec_t &pot) const
{
    std::auto_ptr<Vec_t>	X = checkoutVec();
    std::auto_ptr<Vec_t>	D = checkoutVec();
    std::auto_ptr<Vec_t>	P = checkoutVec();

    X->replicate(src);
    D->replicate(den);
    P->replicate(pot);

    // shuffle
    ShufflePoints(src, *X);
    ShufflePoints(den, *D);
    P->setPointOrder(PointMajor);

    // far interactions
    Error_t status;
    CHK( status = interaction_(*X, *D, *P));

    // shuffle back
    ShufflePoints(*P, pot);

    X->setPointOrder(AxisMajor);	/* ignoring current content */
    D->setPointOrder(AxisMajor);	/* ignoring current content */
    P->setPointOrder(AxisMajor);	/* ignoring current content */

    recycle(X);
    recycle(D);
    recycle(P);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
EvaluateFarInteraction(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const
{
    if ( params_.upsample_interaction ){
	return EvalFarInter_ImpUpsample(src, fi, vel);
    } else {
	return EvalFarInter_Imp(src, fi, vel);
    }
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
EvalFarInter_Imp(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const
{
    std::auto_ptr<Vec_t> den    = checkoutVec();
    std::auto_ptr<Vec_t> slf    = checkoutVec();

    den->replicate(src);
    slf->replicate(vel);

    // multiply area elment into density
    xv(S_.getAreaElement(), fi, *den);

    // compute self (to subtract)
    slf->getDevice().DirectStokes(src.begin(), den->begin(), quad_weights_.begin(),
	slf->getStride(), slf->getNumSubs(), src.begin() /* target */,
	0, slf->getStride() /* number of trgs per surface */,
	slf->begin());

    // incorporating the quadrature weights into the density (use pos as temp holder)
    ax<Sca_t>(quad_weights_, *den, *den);

    CHK(CallInteraction(src, *den, vel));

    // subtract self
    axpy(static_cast<value_type>(-1.0), *slf, vel, vel);

    recycle(den);
    recycle(slf);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
EvalFarInter_ImpUpsample(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const
{
    std::auto_ptr<Vec_t> pos = checkoutVec();
    std::auto_ptr<Vec_t> den = checkoutVec();
    std::auto_ptr<Vec_t> pot = checkoutVec();
    std::auto_ptr<Vec_t> slf = checkoutVec();
    std::auto_ptr<Vec_t> shc = checkoutVec();
    std::auto_ptr<Vec_t> wrk = checkoutVec();

    // prepare for upsampling
    int usf(sht_upsample_.getShOrder());
    pos->resize(src.getNumSubs(), usf);
    den->resize(src.getNumSubs(), usf);
    slf->resize(src.getNumSubs(), usf);
    shc->resize(src.getNumSubs(), usf);
    wrk->resize(src.getNumSubs(), usf);

    // multiply area elment into density (using pot as temp)
    pot->replicate(fi);
    xv(S_.getAreaElement(), fi, *pot);

    // upsample position and density
    Resample( src, sht_, sht_upsample_, *shc, *wrk, *pos);
    Resample(*pot, sht_, sht_upsample_, *shc, *wrk, *den);
    pot->resize(src.getNumSubs(), usf);

    // compute self (to subtract)
    slf->getDevice().DirectStokes(pos->begin(), den->begin(), quad_weights_up_.begin(),
	slf->getStride(), slf->getNumSubs(), pos->begin() /* target */,
	0, slf->getStride() /* number of trgs per surface */,
	slf->begin());

    // incorporating the quadrature weights into the density
    ax<Sca_t>(quad_weights_up_, *den, *den);

    CHK(CallInteraction(*pos, *den, *pot));

    // subtract self
    axpy(static_cast<value_type>(-1.0), *slf, *pot, *slf);

    // downsample
    Resample(*slf, sht_upsample_, sht_, *shc, *wrk, *pot);
    sht_.lowPassFilter(*pot, *wrk, *shc, vel);

    recycle(pos);
    recycle(den);
    recycle(pot);
    recycle(slf);
    recycle(shc);
    recycle(wrk);

    return ErrorEvent::Success;
}

// Linear solve to compute tension such that surface divergence:
// surf_div( velocity + stokes(tension) ) = 0
template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::getTension(
    const Vec_t &vel_in, Sca_t &tension) const
{
    std::auto_ptr<Sca_t> rhs = checkoutSca();
    std::auto_ptr<Sca_t> wrk = checkoutSca();

    S_.div(vel_in, *rhs);

    //! this just negates rhs (not a bug; bad naming for overleaded axpy)
    axpy(static_cast<value_type>(-1), *rhs, *rhs);

    int iter(params_.tension_solver_iter);
    int rsrt(params_.tension_solver_restart);
    value_type tol(params_.tension_solver_tol),relres(params_.tension_solver_tol);
    enum BiCGSReturn solver_ret;
    Error_t ret_val(ErrorEvent::Success);

    COUTDEBUG("Solving for tension");
    solver_ret = linear_solver_(*this, tension, *rhs, rsrt, iter, relres);

    if ( solver_ret  != BiCGSSuccess )
        ret_val = ErrorEvent::DivergenceError;

    COUTDEBUG("Tension solve: Total iter = "<< iter<<", relres = "<<relres);
    COUTDEBUG("Checking true relres");
    ASSERT(((*this)(tension, *wrk),
            axpy(static_cast<value_type>(-1), *wrk, *rhs, *wrk),
            relres = sqrt(AlgebraicDot(*wrk, *wrk))/sqrt(AlgebraicDot(*rhs,*rhs)),
            relres<tol
	    ),
           "relres ("<<relres<<")<tol("<<tol<<")"
           );

    recycle(wrk);
    recycle(rhs);

    return ret_val;
}

// Computes near (self) velocity due to force.
// Computes singular integration on the vesicle surface.
template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::stokes(
    const Vec_t &force, Vec_t &velocity) const
{
    PROFILESTART();

    int imax(S_.getPosition().getGridDim().first);
    int jmax(S_.getPosition().getGridDim().second);
    int np = S_.getPosition().getStride();
    int nv = S_.getPosition().getNumSubs();

    std::auto_ptr<Sca_t> t1 = checkoutSca();
    std::auto_ptr<Sca_t> t2 = checkoutSca();
    std::auto_ptr<Vec_t> v1 = checkoutVec();
    std::auto_ptr<Vec_t> v2 = checkoutVec();

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
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::stokes_double_layer(
    const Vec_t &force, Vec_t &velocity) const
{
    PROFILESTART();

    int imax(S_.getPosition().getGridDim().first);
    int jmax(S_.getPosition().getGridDim().second);
    int np = S_.getPosition().getStride();
    int nv = S_.getPosition().getNumSubs();

    std::auto_ptr<Sca_t> t1 = checkoutSca();
    std::auto_ptr<Sca_t> t2 = checkoutSca();
    std::auto_ptr<Vec_t> v1 = checkoutVec();
    std::auto_ptr<Vec_t> v2 = checkoutVec();
    std::auto_ptr<Vec_t> v3 = checkoutVec();

    ax(w_sph_inv_, S_.getAreaElement(), *t1);

    int numinputs = 4;
    const Sca_t* inputs[] = {&S_.getPosition(), &S_.getNormal(),   &force, t1.get()};
    Sca_t*      outputs[] = { v1.get()        ,  v3.get()      , v2.get(), t2.get()};
    move_pole.setOperands(inputs, numinputs, params_.singular_stokes);

    for(int ii=0;ii < imax; ++ii)
        for(int jj=0;jj < jmax; ++jj)
        {
            move_pole(ii, jj, outputs);

            ax(w_sph_, *t2, *t2);
            xv(*t2, *v2, *v2);

            PROFILESTART();
            S_.getPosition().getDevice().DirectStokesDoubleLayer(v1->begin(), v3->begin(), v2->begin(),
                sing_quad_weights_.begin(), np, nv, S_.getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, velocity.begin());
            PROFILEEND("DblLayerSelfInteraction_",0);
        }

    recycle(t1);
    recycle(t2);
    recycle(v1);
    recycle(v2);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const
{
    std::auto_ptr<Vec_t> fb = checkoutVec();

    COUTDEBUG("Time matvec");
    Intfcl_force_.linearBendingForce(S_, x_new, *fb);
    CHK(stokes(*fb, time_mat_vec));
    axpy(-dt_, time_mat_vec, x_new, time_mat_vec);
    recycle(fb);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::operator()(
    const Sca_t &tension, Sca_t &div_stokes_fs) const
{
    std::auto_ptr<Vec_t> fs = checkoutVec();
    std::auto_ptr<Vec_t> u = checkoutVec();

    COUTDEBUG("Tension matvec");
    Intfcl_force_.tensileForce(S_, tension, *fs);
    CHK(stokes(*fs, *u));
    S_.div(*u, div_stokes_fs);

    recycle(fs);
    recycle(u);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::reparam()
{
    value_type ts(params_.rep_ts);
    value_type vel;

    int ii(-1);
    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();
    std::auto_ptr<Sca_t> wrk = checkoutSca();

    COUTDEBUG("Reparametrization");
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

        COUTDEBUG("Iteration = "<<ii<<", |vel| = "<<vel);

        if(vel < params_.rep_tol )
            break;

    }
    INFO("Iterations = "<<ii<<", |vel| = "<<vel);

    recycle(u1);
    recycle(u2);
    recycle(wrk);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
std::auto_ptr<typename SurfContainer::Sca_t> InterfacialVelocity<SurfContainer, Interaction>::
checkoutSca() const
{
    std::auto_ptr<Sca_t> scp;

    if(scalar_work_q_.empty())
        scp = static_cast<std::auto_ptr<Sca_t> >(new Sca_t);
    else
    {
        scp = static_cast<std::auto_ptr<Sca_t> >(scalar_work_q_.front());
        scalar_work_q_.pop();
    }

    scp->replicate(S_.getPosition());
    ++checked_out_work_sca_;
    return(scp);
}
template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::
recycle(std::auto_ptr<Sca_t> scp) const
{
    scalar_work_q_.push(scp.release());
    --checked_out_work_sca_;
}

template<typename SurfContainer, typename Interaction>
std::auto_ptr<typename SurfContainer::Vec_t> InterfacialVelocity<SurfContainer, Interaction>::
checkoutVec() const
{
    std::auto_ptr<Vec_t> vcp;

    if(vector_work_q_.empty())
        vcp = static_cast<std::auto_ptr<Vec_t> >(new Vec_t);
    else
    {
        vcp = static_cast<std::auto_ptr<Vec_t> >(vector_work_q_.front());
        vector_work_q_.pop();
    }

    vcp->replicate(S_.getPosition());
    ++checked_out_work_vec_;

    return(vcp);
}

template<typename SurfContainer, typename Interaction>
void InterfacialVelocity<SurfContainer, Interaction>::
recycle(std::auto_ptr<Vec_t> vcp) const
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

//////////////////////////////////////////////////////////////////////////
/// DEBUG mode methods ///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
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

    ASSERT(res,"Explicit stepper benchmark failed");
    COUT(emph<<"Explicit stepper benchmark with "
        << S_.getPosition().getNumSubs() << " surface(s) "
        << "Passed\n"
	<< "------------------------------------------------------------"<<emph);
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

    ASSERT(res,"Implicit stepper benchmark failed");
    COUT(emph<<"Implicit stepper benchmark with "
        << S_.getPosition().getNumSubs() << " surface(s) "
        << "Passed\n"
	<< "------------------------------------------------------------"<<emph);

    return ( res );
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkBendingForce(const Vec_t &x, Vec_t &Fb, value_type tol) const
{
    COUTDEBUG("Start benchmarking bending force");
    COUTDEBUG("--------------------------------");

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    Intfcl_force_.linearBendingForce(S_, x, *u1);
    axpy(static_cast<value_type>(-1), Fb, *u1, Fb);

    value_type err = MaxAbs(Fb);
    axpy(static_cast<value_type>(1), *u1, Fb);

    ASSERT(err<tol, "Bending force benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Bending force benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    recycle(u1);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkStokes(const Vec_t &F, Vec_t &SF, value_type tol) const
{
    COUTDEBUG("Start benchmarking singular stokes");
    COUTDEBUG("----------------------------------");

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    stokes(F, *u1);
    axpy(static_cast<value_type>(-1), SF, *u1, SF);

    value_type err = MaxAbs(SF);
    axpy(static_cast<value_type>(1), *u1, SF);

    ASSERT(err<tol, "Singular stokes benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Singular stokes benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    recycle(u1);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkBgFlow(const Vec_t &SFb, Vec_t &vel, value_type tol) const
{
    COUTDEBUG("Start benchmarking background flow");
    COUTDEBUG("----------------------------------");

    //    if (this->interaction_.HasInteraction()){
	this->updateFarField();
	axpy(static_cast<value_type>(1), velocity_, SFb, velocity_);
    // } else {
    // 	axpy(static_cast<value_type>(1), SFb, velocity_);
    // }

    axpy(static_cast<value_type>(-1), vel, velocity_, vel);

    value_type err = MaxAbs(vel);
    axpy(static_cast<value_type>(1), velocity_, vel);

    ASSERT(err<tol, "Background flow benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Background benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkTension(const Vec_t &vel, Sca_t &tension, value_type tol) const
{
    COUTDEBUG("Start benchmarking tension");
    COUTDEBUG("--------------------------");

    std::auto_ptr<Sca_t> wrk = checkoutSca();
    axpy(static_cast<value_type>(0), *wrk, *wrk);
    getTension(vel, *wrk);

    axpy(static_cast<value_type>(-1), tension, *wrk, tension);

    value_type err = MaxAbs(tension);
    axpy(static_cast<value_type>(1), *wrk, tension);

    ASSERT(err<tol, "Tension benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Tension benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    recycle(wrk);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkNewPostitionExplicit(Vec_t &xnew, value_type tol)
{
    COUTDEBUG("Start benchmarking explicit update");
    COUTDEBUG("----------------------------------");

    updateJacobiExplicit(dt_);

    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    value_type err = MaxAbs(xnew);
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);

    ASSERT(err<tol, "Explicit update benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Explicit update benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkTensionImplicit(Sca_t &tension, value_type tol)
{
    COUTDEBUG("Start benchmarking tension (implicit)");
    COUTDEBUG("-------------------------------------");

    this->interaction_.HasInteraction() && this->updateFarField();

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();
    std::auto_ptr<Sca_t> scp = checkoutSca();

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

    ASSERT(err<tol, "Tension (implicit) benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Tension (implicit) benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    recycle(u1);
    recycle(u2);
    recycle(scp);

    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkMatVecImplicit(const Vec_t &x, Vec_t &matvec, value_type tol)
{
    COUTDEBUG("Start benchmarking implicit matvec");
    COUTDEBUG("----------------------------------");

    std::auto_ptr<Vec_t> b = checkoutVec();
    this->operator()(x, *b);

    axpy(static_cast<value_type>(-1), matvec, *b, matvec);
    value_type err = MaxAbs(matvec);
    axpy(static_cast<value_type>(1), *b, matvec);

    ASSERT(err<tol, "Implicit matvec benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Implicit matvec benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    recycle(b);
    return (err < tol);
}

template<typename SurfContainer, typename Interaction>
bool InterfacialVelocity<SurfContainer, Interaction>::
benchmarkNewPostitionImplicit(Vec_t &xnew, value_type tol)
{
    COUTDEBUG("Start benchmarking implicit update");
    COUTDEBUG("----------------------------------");

    updateJacobiGaussSeidel(dt_);

    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    value_type err = MaxAbs(xnew);
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);

    ASSERT(err<tol, "Implicit update benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Implicit update benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    return (err < tol);
}

#endif //NDEBUG
