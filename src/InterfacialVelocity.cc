template<typename SurfContainer, typename Interaction>
InterfacialVelocity<SurfContainer, Interaction>::
InterfacialVelocity(SurfContainer &S_in, const Interaction &Inter,
    const OperatorsMats<Arr_t> &mats,
    const Parameters<value_type> &params, const BgFlowBase<Vec_t> &bgFlow,
    PSolver_t *parallel_solver) :
    S_(S_in),
    interaction_(Inter),
    bg_flow_(bgFlow),
    Intfcl_force_(params,mats),
    params_(params),
    //
    parallel_solver_(parallel_solver),
    psolver_configured_(false),
    precond_configured_(false),
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
    stokes_(mats,params_),
    S_up_(NULL)
{
    pos_vel_.replicate(S_.getPosition());
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

    //spectrum in harmonic space, diagonal
    if (params_.time_precond == DiagonalSpectral){
	position_precond.resize(1,p);
	tension_precond.resize(1,p);
    }
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

    if(S_up_) delete S_up_;
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
updateJacobiExplicit(const SurfContainer& S_, const value_type &dt, Vec_t& dx)
{
    this->dt_ = dt;

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();

    // puts u_inf and interaction in pos_vel_
    this->updateFarField();

    // add S[f_b]
    Intfcl_force_.bendingForce(S_, *u1);
    CHK(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1.0), *u2, pos_vel_, pos_vel_);

    // compute tension
    CHK(getTension(pos_vel_, tension_));

    // add S[f_sigma]
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    CHK(stokes(*u1, *u2));
    axpy(static_cast<value_type>(1.0), *u2, pos_vel_, pos_vel_);

    //axpy(dt_, pos_vel_, S_.getPosition(), S_.getPositionModifiable());
    dx.replicate(S_.getPosition());
    axpy(dt, pos_vel_, dx);

    recycle(u1);
    recycle(u2);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateJacobiGaussSeidel(const SurfContainer& S_, const value_type &dt, Vec_t& dx)
{
    this->dt_ = dt;

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();
    std::auto_ptr<Vec_t> u3 = checkoutVec();

    // put far field in pos_vel_ and the sum with S[f_b] in u1
    this->updateFarField();
    Intfcl_force_.bendingForce(S_, *u2);
    CHK(stokes(*u2, *u1));
    axpy(static_cast<value_type>(1.0), pos_vel_, *u1, *u1);

    // tension
    CHK(getTension(*u1, tension_));
    Intfcl_force_.tensileForce(S_, tension_, *u1);
    CHK(stokes(*u1, *u2));

    // position rhs
    axpy(static_cast<value_type>(1.0), pos_vel_, *u2, *u1);
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

    //u2->getDevice().Memcpy(S_.getPositionModifiable().begin(), u2->begin(),
    //    S_.getPosition().size() * sizeof(value_type),
    //    u2->getDevice().MemcpyDeviceToDevice);
    dx.replicate(S_.getPosition());
    axpy(-1, S_.getPosition(), *u2, dx);

    recycle(u1);
    recycle(u2);
    recycle(u3);

    return ret_val;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateImplicit(const SurfContainer& S_, const value_type &dt, Vec_t& dx)
{
    PROFILESTART();
    this->dt_ = dt;
    SolverScheme scheme(GloballyImplicit);
    INFO("Taking a time step using "<<scheme<<" scheme");
    CHK(Prepare(scheme));

    if (params_.solve_for_velocity) {
	CHK(AssembleRhsVel(parallel_rhs_, dt_, scheme));
    } else {
	CHK(AssembleRhsPos(parallel_rhs_, dt_, scheme));
    }

    Error_t err=ErrorEvent::Success;
    if(err==ErrorEvent::Success) err=AssembleInitial(parallel_u_, dt_, scheme);
    if(err==ErrorEvent::Success) err=Solve(parallel_rhs_, parallel_u_, dt_, scheme);
    if(err==ErrorEvent::Success) err=Update(parallel_u_);

    dx.replicate(S_.getPosition());
    if (params_.solve_for_velocity){
	//axpy(dt_, pos_vel_, S_.getPosition(), S_.getPositionModifiable());
        axpy(dt, pos_vel_, dx);
    } else {
	//S_.setPosition(pos_vel_);
        axpy(-1, S_.getPosition(), pos_vel_, dx);
    }

    PROFILEEND("",0);
    return err;
}

template<typename SurfContainer, typename Interaction>
size_t InterfacialVelocity<SurfContainer, Interaction>::stokesBlockSize() const{

    return (params_.pseudospectral ?
	S_.getPosition().size() :
	S_.getPosition().getShOrder()*(S_.getPosition().getShOrder()+2)*S_.getPosition().getNumSubFuncs() ); /* (p+1)^2-1 (last freq doesn't have a cosine */
}

template<typename SurfContainer, typename Interaction>
size_t InterfacialVelocity<SurfContainer, Interaction>::tensionBlockSize() const{
    return stokesBlockSize()/3;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::Prepare(const SolverScheme &scheme) const
{
    PROFILESTART();

    if (pos_vel_.size() != S_.getPosition().size() ||
	dl_coeff_.size() != S_.getPosition().getNumSubs() ){

      COUTDEBUG("Resizing the containers");
      pos_vel_.replicate(S_.getPosition());
      tension_.replicate(S_.getPosition());

      COUTDEBUG("zeroing content of velocity and tension arrays");
      pos_vel_.getDevice().Memset(pos_vel_.begin(), 0, sizeof(value_type)*pos_vel_.size());
      tension_.getDevice().Memset(tension_.begin(), 0, sizeof(value_type)*tension_.size());

      //Permitting the viscosity constrast to be different for each vesicle
      size_t     nves(pos_vel_.getNumSubs());
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

    ASSERT(   pos_vel_.size() ==   S_.getPosition().size(), "inccorrect size");
    ASSERT( 3*tension_.size() ==   S_.getPosition().size(), "inccorrect size");
    ASSERT(  dl_coeff_.size() ==   S_.getPosition().getNumSubs(), "inccorrect size");
    ASSERT( vel_coeff_.size() ==   S_.getPosition().getNumSubs(), "inccorrect size");

    INFO("Setting interaction source and target");
    stokes_.SetSrcCoord(S_);
    stokes_.SetTrgCoord(S_);

    if (!precond_configured_ && params_.time_precond!=NoPrecond)
	ConfigurePrecond(params_.time_precond);

    //!@bug doesn't support repartitioning
    if (!psolver_configured_ && scheme==GloballyImplicit){
	ASSERT(parallel_solver_ != NULL, "need a working parallel solver");
	CHK(ConfigureSolver(scheme));
    }

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ConfigureSolver(const SolverScheme &scheme) const
{
    PROFILESTART();
    ASSERT(scheme==GloballyImplicit, "Unsupported scheme");
    COUTDEBUG("Configuring the parallel solver");

    typedef typename PSolver_t::matvec_type POp;
    typedef typename PSolver_t::vec_type PVec;
    typedef typename PVec::size_type size_type;

    // Setting up the operator
    size_type sz(stokesBlockSize() + tensionBlockSize());
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
    CHK(parallel_solver_->SetTolerances(params_.time_tol,
	    PSolver_t::PLS_DEFAULT,
	    PSolver_t::PLS_DEFAULT,
	    params_.time_iter_max));

    CHK(parallel_solver_->Configure());

    // setting up the preconditioner
    if (params_.time_precond != NoPrecond){
	ASSERT(precond_configured_, "The preconditioner isn't configured yet");
	CHK(parallel_solver_->SetPrecondContext(static_cast<const void*>(this)));
	CHK(parallel_solver_->UpdatePrecond(ImplicitPrecond));
    }
    psolver_configured_ = true;

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ConfigurePrecond(const PrecondScheme &precond) const{

    PROFILESTART();
    if (precond!=DiagonalSpectral)
	return ErrorEvent::NotImplementedError; /* Unsupported preconditioner scheme */

    INFO("Setting up the diagonal preceonditioner");
    value_type *buffer = new value_type[position_precond.size() * sizeof(value_type)];

    { //bending precond
	int idx(0), N(0);
	// The sh coefficients are ordered by m and then n
	for(int iM=0; iM<position_precond.getGridDim().second; ++iM){
	    for(int iN=++N/2; iN<position_precond.getGridDim().first; ++iN){
		value_type bending_precond(1.0/fabs(1.0-dt_*iN*iN*iN));
		bending_precond = fabs(bending_precond) < 1e3   ? bending_precond : 1.0;
		buffer[idx]     = fabs(bending_precond) > 1e-10 ? bending_precond : 1.0;
		++idx;
	    }
	}
	position_precond.getDevice().Memcpy(position_precond.begin(), buffer,
	    position_precond.size() * sizeof(value_type),
	    device_type::MemcpyHostToDevice);
    }
    { // tension precond
	int idx(0), N(0);
	for(int iM=0; iM<tension_precond.getGridDim().second; ++iM){
	    for(int iN=++N/2; iN<tension_precond.getGridDim().first; ++iN){
		value_type eig(4*iN*iN-1);
		eig *= 2*iN+3;
		eig /= iN+1;
		eig /= 2*iN*iN+2*iN-1;
		eig  = iN==0 ? 1.0 : eig/iN;
		buffer[idx] = fabs(eig) > 1e-10 ? eig : 1.0;
		++idx;
	    }
	}
	tension_precond.getDevice().Memcpy(tension_precond.begin(), buffer,
	    tension_precond.size() * sizeof(value_type),
	    device_type::MemcpyHostToDevice);
    }

    delete[] buffer;
    precond_configured_=true;
    PROFILEEND("",0);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
AssembleRhsVel(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const
{
    PROFILESTART();
    ASSERT(scheme==GloballyImplicit, "Unsupported scheme");
    INFO("Assembling RHS to solve for velocity");

    // rhs=[u_inf+Bx;div(u_inf+Bx)]
    COUTDEBUG("Evaluate background flow");
    std::auto_ptr<Vec_t> vRhs = checkoutVec();
    vRhs->replicate(S_.getPosition());
    CHK(BgFlow(*vRhs, dt));

    COUTDEBUG("Computing the far-field interaction due to explicit traction jump");
    std::auto_ptr<Vec_t> f  = checkoutVec();
    std::auto_ptr<Vec_t> Sf = checkoutVec();
    Intfcl_force_.explicitTractionJump(S_, *f);
    stokes_.SetDensitySL(f.get());
    stokes_.SetDensityDL(NULL);
    stokes_(*Sf);
    axpy(static_cast<value_type>(1.0), *Sf, *vRhs, *vRhs);

    COUTDEBUG("Computing rhs for div(u)");
    std::auto_ptr<Sca_t> tRhs = checkoutSca();
    S_.div(*vRhs, *tRhs);

    ASSERT( vRhs->getDevice().isNumeric(vRhs->begin(), vRhs->size()), "Non-numeric rhs");
    ASSERT( tRhs->getDevice().isNumeric(tRhs->begin(), tRhs->size()), "Non-numeric rhs");

    // copy data
    size_t xsz(stokesBlockSize()), tsz(tensionBlockSize());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;
    CHK(parallel_rhs_->GetArray(i, rsz));
    ASSERT(rsz==xsz+tsz,"Bad sizes");

    if (params_.pseudospectral){
	COUTDEBUG("Copy data to parallel rhs array");
	vRhs->getDevice().Memcpy(i    , vRhs->begin(), xsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tRhs->getDevice().Memcpy(i+xsz, tRhs->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    } else {  /* Galerkin */
	COUTDEBUG("Project RHS to spectral coefficient");
	std::auto_ptr<Vec_t> vRhsSh  = checkoutVec();
	std::auto_ptr<Sca_t> tRhsSh  = checkoutSca();
	std::auto_ptr<Vec_t> wrk     = checkoutVec();

	vRhsSh->replicate(*vRhs);
	tRhsSh->replicate(*tRhs);
	wrk->replicate(*vRhs);

	sht_.forward(*vRhs, *wrk, *vRhsSh);
	sht_.forward(*tRhs, *wrk, *tRhsSh);

	COUTDEBUG("Copy data to parallel RHS array (size="<<xsz<<"+"<<tsz<<")");
	vRhs->getDevice().Memcpy(i    , vRhsSh->begin(), xsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tRhs->getDevice().Memcpy(i+xsz, tRhsSh->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);

	recycle(vRhsSh);
	recycle(tRhsSh);
	recycle(wrk);
    }

    CHK(parallel_rhs_->RestoreArray(i));

    recycle(vRhs);
    recycle(f);
    recycle(Sf);
    recycle(tRhs);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
AssembleRhsPos(PVec_t *rhs, const value_type &dt, const SolverScheme &scheme) const
{
    PROFILESTART();
    ASSERT(scheme==GloballyImplicit, "Unsupported scheme");
    INFO("Assembling RHS to solve for position");

    COUTDEBUG("Evaluate background flow");
    std::auto_ptr<Vec_t> pRhs = checkoutVec();
    std::auto_ptr<Vec_t> pRhs2 = checkoutVec();
    pRhs->replicate(S_.getPosition());
    pRhs2->replicate(S_.getPosition());
    CHK(BgFlow(*pRhs, dt));

    if( fabs(params_.viscosity_contrast-1.0)>1e-12){
	COUTDEBUG("Computing the rhs due to viscosity difference");
	std::auto_ptr<Vec_t> x  = checkoutVec();
	std::auto_ptr<Vec_t> Dx = checkoutVec();
	av(dl_coeff_, S_.getPosition(), *x);
	stokes_.SetDensitySL(NULL);
	stokes_.SetDensityDL(x.get());
	stokes_(*Dx);
	axpy(-dt, *pRhs, *Dx, *pRhs);
	axpy(static_cast<value_type>(-1.0), *pRhs, *pRhs);

	recycle(x);
	recycle(Dx);
    } else
	axpy(dt, *pRhs, *pRhs);

    COUTDEBUG("Computing rhs for div(u)");
    std::auto_ptr<Sca_t> tRhs = checkoutSca();
    S_.div(*pRhs, *tRhs);

    av(vel_coeff_, S_.getPosition(), *pRhs2);
    axpy(static_cast<value_type>(1.0), *pRhs, *pRhs2, *pRhs);

    ASSERT( pRhs->getDevice().isNumeric(pRhs->begin(), pRhs->size()), "Non-numeric rhs");
    ASSERT( tRhs->getDevice().isNumeric(tRhs->begin(), tRhs->size()), "Non-numeric rhs");

    // copy data
    size_t xsz(stokesBlockSize()), tsz(tensionBlockSize());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;
    CHK(parallel_rhs_->GetArray(i, rsz));
    ASSERT(rsz==xsz+tsz,"Bad sizes");

    if (params_.pseudospectral){
	COUTDEBUG("Copy data to parallel rhs array");
	pRhs->getDevice().Memcpy(i    , pRhs->begin(), xsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tRhs->getDevice().Memcpy(i+xsz, tRhs->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    } else {  /* Galerkin */
	COUTDEBUG("Project RHS to spectral coefficient");
	std::auto_ptr<Vec_t> pRhsSh  = checkoutVec();
	std::auto_ptr<Sca_t> tRhsSh  = checkoutSca();
	std::auto_ptr<Vec_t> wrk     = checkoutVec();

	pRhsSh->replicate(*pRhs);
	tRhsSh->replicate(*tRhs);
	wrk->replicate(*pRhs);

	sht_.forward(*pRhs, *wrk, *pRhsSh);
	sht_.forward(*tRhs, *wrk, *tRhsSh);

	COUTDEBUG("Copy data to parallel rhs array");
	pRhs->getDevice().Memcpy(i    , pRhsSh->begin(), xsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tRhs->getDevice().Memcpy(i+xsz, tRhsSh->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);

	recycle(pRhsSh);
	recycle(tRhsSh);
	recycle(wrk);
    }

    CHK(parallel_rhs_->RestoreArray(i));

    recycle(pRhs);
    recycle(pRhs2);
    recycle(tRhs);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
AssembleInitial(PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const
{
    PROFILESTART();
    COUTDEBUG("Using current position/tension as initial guess");
    size_t vsz(stokesBlockSize()), tsz(tensionBlockSize());
    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;

    ASSERT( pos_vel_.getDevice().isNumeric(pos_vel_.begin(), pos_vel_.size()), "Non-numeric velocity");
    ASSERT( tension_.getDevice().isNumeric(tension_.begin(), tension_.size()), "Non-numeric tension");

    CHK(parallel_u_->GetArray(i, rsz));
    ASSERT(rsz==vsz+tsz,"Bad sizes");

    if (params_.pseudospectral){
	COUTDEBUG("Copy initial guess to parallel solution array");
	pos_vel_.getDevice().Memcpy(i    , pos_vel_.begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tension_.getDevice().Memcpy(i+vsz, tension_.begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    } else {  /* Galerkin */
    	COUTDEBUG("Project initial guess to spectral coefficient");
	std::auto_ptr<Vec_t> voxSh  = checkoutVec();
	std::auto_ptr<Sca_t> tSh    = checkoutSca();
	std::auto_ptr<Vec_t> wrk    = checkoutVec();

	voxSh->replicate(pos_vel_);
	tSh->replicate(tension_);
	wrk->replicate(pos_vel_);

	sht_.forward(pos_vel_, *wrk, *voxSh);
	sht_.forward(tension_, *wrk, *tSh);

	COUTDEBUG("Copy initial guess to parallel solution array");
	voxSh->getDevice().Memcpy(i    , voxSh->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tSh->getDevice().Memcpy(  i+vsz, tSh->begin()  , tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);

	recycle(voxSh);
	recycle(tSh);
	recycle(wrk);
    }

    CHK(parallel_u_->RestoreArray(i));
    CHK(parallel_solver_->InitialGuessNonzero(true));
    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::ImplicitMatvecPhysical(Vec_t &vox, Sca_t &ten) const
{
    PROFILESTART();

    std::auto_ptr<Vec_t> f   = checkoutVec();
    std::auto_ptr<Vec_t> Sf  = checkoutVec();
    std::auto_ptr<Vec_t> Du  = checkoutVec();
    f->replicate(vox);
    Sf->replicate(vox);
    Du->replicate(vox);

    COUTDEBUG("Computing the interfacial forces and setting single-layer density");
    if (params_.solve_for_velocity){
	// Bending of dt*u + tension of sigma
	axpy(dt_, vox, *Du);
	Intfcl_force_.implicitTractionJump(S_, *Du, ten, *f);
    } else {
	// dt*(Bending of x + tension of sigma)
	Intfcl_force_.implicitTractionJump(S_, vox, ten, *f);
	axpy(dt_, *f, *f);
    }
    stokes_.SetDensitySL(f.get());

    if( fabs(params_.viscosity_contrast-1.0)>1e-12){
	COUTDEBUG("Setting the double-layer density");
	av(dl_coeff_, vox, *Du);
	stokes_.SetDensityDL(Du.get());
    } else {
	stokes_.SetDensityDL(NULL);
    }

    COUTDEBUG("Calling stokes");
    stokes_(*Sf);

    COUTDEBUG("Computing the div term");
    //! @note For some reason, doing the linear algebraic manipulation
    //! and writing the constraint as -\div{S[f_b+f_s]} = \div{u_inf
    //! almost halves the number of gmres iterations. Also having the
    //! minus sign in the matvec is tangibly better (1-2
    //! iterations). Need to investigate why.
    S_.div(*Sf, ten);
    axpy((value_type) -1.0, ten, ten);

    if( fabs(params_.viscosity_contrast-1.0)>1e-12)
	av(vel_coeff_, vox, vox);

    axpy((value_type) -1.0, *Sf, vox, vox);

    ASSERT(vox.getDevice().isNumeric(vox.begin(), vox.size()), "Non-numeric velocity");
    ASSERT(ten.getDevice().isNumeric(ten.begin(), ten.size()), "Non-numeric divergence");

    recycle(f);
    recycle(Sf);
    recycle(Du);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ImplicitApply(const POp_t *o, const value_type *x, value_type *y)
{
    PROFILESTART();
    const InterfacialVelocity *F(NULL);
    o->Context((const void**) &F);
    size_t vsz(F->stokesBlockSize()), tsz(F->tensionBlockSize());

    std::auto_ptr<Vec_t> vox = F->checkoutVec();
    std::auto_ptr<Sca_t> ten = F->checkoutSca();
    vox->replicate(F->pos_vel_);
    ten->replicate(F->tension_);

    COUTDEBUG("Unpacking the input from parallel vector");
    if (F->params_.pseudospectral){
	vox->getDevice().Memcpy(vox->begin(), x    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	ten->getDevice().Memcpy(ten->begin(), x+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    } else {  /* Galerkin */
	std::auto_ptr<Vec_t> voxSh = F->checkoutVec();
	std::auto_ptr<Sca_t> tSh   = F->checkoutSca();
	std::auto_ptr<Vec_t> wrk   = F->checkoutVec();

	voxSh->replicate(*vox);
	tSh->replicate(*ten);
	wrk->replicate(*vox);
	voxSh->getDevice().Memcpy(voxSh->begin(), x    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	tSh  ->getDevice().Memcpy(tSh->begin()  , x+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);

	COUTDEBUG("Mapping the input to physical space");
	F->sht_.backward(*voxSh, *wrk, *vox);
	F->sht_.backward(*tSh  , *wrk, *ten);

	F->recycle(voxSh);
	F->recycle(tSh);
	F->recycle(wrk);
    }

    F->ImplicitMatvecPhysical(*vox, *ten);

    if (F->params_.pseudospectral){
	COUTDEBUG("Packing the matvec into parallel vector");
	vox->getDevice().Memcpy(y    , vox->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	ten->getDevice().Memcpy(y+vsz, ten->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    } else {  /* Galerkin */
	COUTDEBUG("Mapping the matvec to physical space");
	std::auto_ptr<Vec_t> voxSh = F->checkoutVec();
	std::auto_ptr<Sca_t> tSh   = F->checkoutSca();
	std::auto_ptr<Vec_t> wrk   = F->checkoutVec();

	voxSh->replicate(*vox);
	tSh->replicate(*ten);
	wrk->replicate(*vox);

	F->sht_.forward(*vox, *wrk, *voxSh);
	F->sht_.forward(*ten, *wrk, *tSh);

	COUTDEBUG("Packing the matvec into parallel vector");
	voxSh->getDevice().Memcpy(y    , voxSh->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tSh  ->getDevice().Memcpy(y+vsz, tSh->begin()  , tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);

	F->recycle(voxSh);
	F->recycle(tSh);
	F->recycle(wrk);
    }

    F->recycle(vox);
    F->recycle(ten);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
ImplicitPrecond(const PSolver_t *ksp, const value_type *x, value_type *y)
{
    PROFILESTART();
    const InterfacialVelocity *F(NULL);
    ksp->PrecondContext((const void**) &F);

    size_t vsz(F->stokesBlockSize()), tsz(F->tensionBlockSize());

    std::auto_ptr<Vec_t> vox = F->checkoutVec();
    std::auto_ptr<Vec_t> vxs = F->checkoutVec();
    std::auto_ptr<Vec_t> wrk = F->checkoutVec();
    vox->replicate(F->pos_vel_);
    vxs->replicate(F->pos_vel_);
    wrk->replicate(F->pos_vel_);

    std::auto_ptr<Sca_t> ten = F->checkoutSca();
    std::auto_ptr<Sca_t> tns = F->checkoutSca();
    ten->replicate(F->tension_);
    tns->replicate(F->tension_);

    COUTDEBUG("Unpacking the input parallel vector");
    if (F->params_.pseudospectral){
	vox->getDevice().Memcpy(vox->begin(), x    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	ten->getDevice().Memcpy(ten->begin(), x+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	F->sht_.forward(*vox, *wrk, *vxs);
	F->sht_.forward(*ten, *wrk, *tns);
    } else {  /* Galerkin */
	vxs->getDevice().Memcpy(vxs->begin(), x    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	tns->getDevice().Memcpy(tns->begin(), x+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    }

    COUTDEBUG("Applying diagonal preconditioner");
    F->sht_.ScaleFreq(vxs->begin(), vxs->getNumSubFuncs(), F->position_precond.begin(), vxs->begin());
    F->sht_.ScaleFreq(tns->begin(), tns->getNumSubFuncs(), F->tension_precond.begin() , tns->begin());

    if (F->params_.pseudospectral){
	F->sht_.backward(*vxs, *wrk, *vox);
	F->sht_.backward(*tns, *wrk, *ten);
	vox->getDevice().Memcpy(y    , vox->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	ten->getDevice().Memcpy(y+vsz, ten->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    } else {  /* Galerkin */
	vxs->getDevice().Memcpy(y    , vxs->begin(), vsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
	tns->getDevice().Memcpy(y+vsz, tns->begin(), tsz * sizeof(value_type), device_type::MemcpyDeviceToHost);
    }

    F->recycle(vox);
    F->recycle(vxs);
    F->recycle(wrk);
    F->recycle(ten);
    F->recycle(tns);

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
Solve(const PVec_t *rhs, PVec_t *u0, const value_type &dt, const SolverScheme &scheme) const
{
    PROFILESTART();
    INFO("Solving for position/velocity and tension using "<<scheme<<" scheme.");

    CHK(parallel_solver_->Solve(parallel_rhs_, parallel_u_));
    typename PVec_t::size_type	iter;
    CHK(parallel_solver_->IterationNumber(iter));
    INFO("Parallel solver returned after "<<iter<<" iteration(s).");
    INFO("Parallal solver report:");
    Error_t err=parallel_solver_->ViewReport();

    PROFILEEND("",0);
    return err;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::Update(PVec_t *u0)
{
    PROFILESTART();
    COUTDEBUG("Updating position and tension.");
    size_t vsz(stokesBlockSize()), tsz(tensionBlockSize());

    typename PVec_t::iterator i(NULL);
    typename PVec_t::size_type rsz;

    CHK(u0->GetArray(i, rsz));
    ASSERT(rsz==vsz+tsz,"Bad sizes");

    if (params_.pseudospectral){
	COUTDEBUG("Copy data from parallel solution array");
	pos_vel_.getDevice().Memcpy(pos_vel_.begin(), i    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	tension_.getDevice().Memcpy(tension_.begin(), i+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);
    } else { /* Galerkin */
	COUTDEBUG("Unpacking the solution from parallel vector");
	std::auto_ptr<Vec_t> voxSh = checkoutVec();
	std::auto_ptr<Sca_t> tSh   = checkoutSca();
	std::auto_ptr<Vec_t> wrk   = checkoutVec();

	voxSh->replicate(pos_vel_);
	tSh->replicate(tension_);
	wrk->replicate(pos_vel_);

	voxSh->getDevice().Memcpy(voxSh->begin(), i    , vsz * sizeof(value_type), device_type::MemcpyHostToDevice);
	tSh  ->getDevice().Memcpy(tSh  ->begin(), i+vsz, tsz * sizeof(value_type), device_type::MemcpyHostToDevice);

	COUTDEBUG("Mapping the solution to physical space");
	sht_.backward(*voxSh, *wrk, pos_vel_);
	sht_.backward(*tSh  , *wrk, tension_);

	recycle(voxSh);
	recycle(tSh);
	recycle(wrk);
    }

    CHK(u0->RestoreArray(i));

    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
BgFlow(Vec_t &bg, const value_type &dt) const{
    //!@bug the time should be passed to the BgFlow handle.
    bg_flow_(S_.getPosition(), 0, bg);

    return ErrorEvent::Success;
}

// Compute velocity_far = velocity_bg + FMM(bending+tension) - DirectStokes(bending+tension)
template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
updateFarField() const
{
    PROFILESTART();
    pos_vel_.replicate(S_.getPosition());
    CHK(this->BgFlow(pos_vel_, this->dt_));

    if (this->interaction_.HasInteraction()){
	std::auto_ptr<Vec_t>	fi  = checkoutVec();
	std::auto_ptr<Vec_t>	vel = checkoutVec();
	fi->replicate(pos_vel_);
	vel->replicate(pos_vel_);

	Intfcl_force_.bendingForce(S_, *fi);
	Intfcl_force_.tensileForce(S_, tension_, *vel);
	axpy(static_cast<value_type>(1.0), *fi, *vel, *fi);

	EvaluateFarInteraction(S_.getPosition(), *fi, *vel);
	axpy(static_cast<value_type>(1.0), *vel, pos_vel_, pos_vel_);

	recycle(fi);
	recycle(vel);
    }
    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
CallInteraction(const Vec_t &src, const Vec_t &den, Vec_t &pot) const
{
    PROFILESTART();
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
    PROFILEEND("",0);

    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
EvaluateFarInteraction(const Vec_t &src, const Vec_t &fi, Vec_t &vel) const
{
    if ( params_.interaction_upsample ){
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
	slf->getStride(), slf->getStride(), slf->getNumSubs(), src.begin() /* target */,
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
	slf->getStride(), slf->getStride(), slf->getNumSubs(), pos->begin() /* target */,
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
    PROFILESTART();
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
    PROFILEEND("",0);

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

            S_.getPosition().getDevice().DirectStokes(v1->begin(), v2->begin(),
                sing_quad_weights_.begin(), np, np, nv, S_.getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, velocity.begin());
        }

    recycle(t1);
    recycle(t2);
    recycle(v1);
    recycle(v2);

    PROFILEEND("SelfInteraction_",0);
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

            S_.getPosition().getDevice().DirectStokesDoubleLayer(v1->begin(), v3->begin(), v2->begin(),
                sing_quad_weights_.begin(), np, nv, S_.getPosition().begin(),
                ii * jmax + jj, ii * jmax + jj + 1, velocity.begin());

        }

    recycle(t1);
    recycle(t2);
    recycle(v1);
    recycle(v2);

    PROFILEEND("DblLayerSelfInteraction_",0);
    return ErrorEvent::Success;
}

template<typename SurfContainer, typename Interaction>
Error_t InterfacialVelocity<SurfContainer, Interaction>::
operator()(const Vec_t &x_new, Vec_t &time_mat_vec) const
{
    PROFILESTART();
    std::auto_ptr<Vec_t> fb = checkoutVec();

    COUTDEBUG("Time matvec");
    Intfcl_force_.linearBendingForce(S_, x_new, *fb);
    CHK(stokes(*fb, time_mat_vec));
    axpy(-dt_, time_mat_vec, x_new, time_mat_vec);
    recycle(fb);
    PROFILEEND("",0);

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
    PROFILESTART();

    value_type ts(params_.rep_ts);
    value_type vel(params_.rep_tol+1);
    value_type last_vel(vel+1);
    int flag(1);

    int ii(-1);
    SurfContainer* Surf;

    std::auto_ptr<Vec_t> u1 = checkoutVec();
    std::auto_ptr<Vec_t> u2 = checkoutVec();
    std::auto_ptr<Sca_t> wrk = checkoutSca();

    if (params_.rep_upsample){
        INFO("Upsampling for reparametrization");
        S_.resample(params_.upsample_freq, &S_up_);
        Surf = S_up_;
    } else {
        Surf = &S_;
    }

    u1 ->replicate(Surf->getPosition());
    u2 ->replicate(Surf->getPosition());
    wrk->replicate(Surf->getPosition());

    std::vector<value_type> pole_quad;
    { // compute quadrature to find pole
      size_t p=Surf->getShOrder();

      //Gauss-Legendre quadrature nodes and weights
      std::vector<double> x(p+1),w(p+1);
      cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

      std::vector<double> leg((p+1)*(p+1));
      LegPoly(&x[0],x.size(),p+1,&leg[0]);
      pole_quad.resize(p+1);
      for(size_t i=0;i<pole_quad.size();i++) pole_quad[i]=0;
      for(size_t j=0;j<p+1;j++){
        for(size_t i=0;i<p+1;i++){
          pole_quad[i]+=leg[j*(p+1)+i]*sqrt(2.0*j+1.0);
        }
      }
      for(size_t i=0;i<p+1;i++){
        pole_quad[i]*=w[i]*0.25/p;
      }
    }

    ii=-1;
    while ( ++ii < params_.rep_maxit )
    { // rearrange mesh points in nearly regular grid.
        { // compute u1  (break loop if mesh is good)
          const Vec_t& v=Surf->getPosition();
          int sh_order=v.getShOrder();
          size_t omp_p=omp_get_max_threads();
          size_t N_ves   =v.getNumSubs(); // Number of vesicles
          size_t M_ves   =(1+sh_order)*(2*sh_order);
          size_t M_ves_  =(3+sh_order)*(2*sh_order);

          std::vector<value_type> vx(N_ves*M_ves_);
          std::vector<value_type> vy(N_ves*M_ves_);
          std::vector<value_type> vz(N_ves*M_ves_);

          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){
            size_t a=((tid+0)*N_ves)/omp_p;
            size_t b=((tid+1)*N_ves)/omp_p;
            for(size_t i=a;i<b;i++){
              const value_type* x=v.getSubN_begin(i)+0*M_ves;
              const value_type* y=v.getSubN_begin(i)+1*M_ves;
              const value_type* z=v.getSubN_begin(i)+2*M_ves;

              for(size_t j=0;j<M_ves;j++){
                vx[i*M_ves_+j +2*sh_order]=x[j];
                vy[i*M_ves_+j +2*sh_order]=y[j];
                vz[i*M_ves_+j +2*sh_order]=z[j];
              }

              value_type pole[2*COORD_DIM];
              for(size_t j=0;j<2*COORD_DIM;j++) pole[j]=0;
              for(size_t k0=0;k0<(1+sh_order);k0++){
                for(size_t k1=0;k1<(2*sh_order);k1++){
                  size_t k=k1+k0*(2*sh_order);
                  pole[0*COORD_DIM+0]+=pole_quad[sh_order-k0]*x[k];
                  pole[0*COORD_DIM+1]+=pole_quad[sh_order-k0]*y[k];
                  pole[0*COORD_DIM+2]+=pole_quad[sh_order-k0]*z[k];
                  pole[1*COORD_DIM+0]+=pole_quad[         k0]*x[k];
                  pole[1*COORD_DIM+1]+=pole_quad[         k0]*y[k];
                  pole[1*COORD_DIM+2]+=pole_quad[         k0]*z[k];
                }
              }
              for(size_t j=0;j<2*sh_order;j++){
                vx[i*M_ves_+       j  ]=pole[0*COORD_DIM+0];
                vy[i*M_ves_+       j  ]=pole[0*COORD_DIM+1];
                vz[i*M_ves_+       j  ]=pole[0*COORD_DIM+2];
                vx[i*M_ves_+M_ves_-j-1]=pole[1*COORD_DIM+0];
                vy[i*M_ves_+M_ves_-j-1]=pole[1*COORD_DIM+1];
                vz[i*M_ves_+M_ves_-j-1]=pole[1*COORD_DIM+2];
              }
            }
            for(size_t i=a;i<b;i++){
              value_type* x=u1->getSubN_begin(i)+0*M_ves;
              value_type* y=u1->getSubN_begin(i)+1*M_ves;
              value_type* z=u1->getSubN_begin(i)+2*M_ves;
              for(size_t k0=0;k0<(1+sh_order);k0++){
                for(size_t k1=0;k1<(2*sh_order);k1++){
                  size_t k=k1+k0*(2*sh_order);
                  x[k]=0;
                  y[k]=0;
                  z[k]=0;

                  size_t k1_0=((2*sh_order+k1-1)%(2*sh_order));
                  size_t k1_1=((2*sh_order+k1+1)%(2*sh_order));



                  value_type dx;
                  value_type dy;
                  value_type dz;

                  dx=(vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dy=(vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dz=(vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r0=sqrt(dx*dx+dy*dy+dz*dz);

                  dx=(vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dy=(vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dz=(vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r1=sqrt(dx*dx+dy*dy+dz*dz);

                  dx=(vx[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dy=(vy[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dz=(vz[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r2=sqrt(dx*dx+dy*dy+dz*dz);

                  dx=(vx[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dy=(vy[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  dz=(vz[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r3=sqrt(dx*dx+dy*dy+dz*dz);



                  value_type r=r0+r1+r2+r3;

                  x[k]+=vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]*r0/r;
                  y[k]+=vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]*r0/r;
                  z[k]+=vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1_0]*r0/r;

                  x[k]+=vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]*r1/r;
                  y[k]+=vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]*r1/r;
                  z[k]+=vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]*r1/r;

                  x[k]+=vx[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]*r2/r;
                  y[k]+=vy[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]*r2/r;
                  z[k]+=vz[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]*r2/r;

                  x[k]+=vx[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]*r3/r;
                  y[k]+=vy[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]*r3/r;
                  z[k]+=vz[i*M_ves_ +(k0+2)*(2*sh_order)+k1  ]*r3/r;
                }
              }
            }
          }

          { // check if mesh is good
            value_type max_ratio=1.0;
            for(size_t i=0;i<N_ves;i++){
              for(size_t k0=0;k0<(1+sh_order);k0++){
                value_type max_=0   ;
                value_type min_=1e10;
                for(size_t k1=0;k1<(2*sh_order);k1++){
                  size_t k=k1+k0*(2*sh_order);
                  size_t k1_1=((2*sh_order+k1+1)%(2*sh_order));
                  value_type dx=(vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type dy=(vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type dz=(vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1_1]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r=sqrt(dx*dx+dy*dy+dz*dz);
                  max_=std::max(max_,r);
                  min_=std::min(min_,r);
                }
                max_ratio=std::max(max_ratio,max_/min_);
              }
              for(size_t k1=0;k1<(2*sh_order);k1++){
                value_type max_=0   ;
                value_type min_=1e10;
                for(size_t k0=0;k0<=(1+sh_order);k0++){
                  size_t k=k1+k0*(2*sh_order);
                  value_type dx=(vx[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vx[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type dy=(vy[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vy[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type dz=(vz[i*M_ves_ +(k0+0)*(2*sh_order)+k1  ]-vz[i*M_ves_ +(k0+1)*(2*sh_order)+k1]);
                  value_type r=sqrt(dx*dx+dy*dy+dz*dz);
                  max_=std::max(max_,r);
                  min_=std::min(min_,r);
                }
                max_ratio=std::max(max_ratio,max_/min_);
              }
            }
            COUTDEBUG("Iteration = "<<ii<<", Mesh-irregularity = "<<max_ratio);
            if(ii==0 && max_ratio<4) break;
            else if(ii>0 && max_ratio<2) break;
          }

          { // project to sph space
            Vec_t wrk[3];
            wrk[0].resize(v.getNumSubs(), sh_order);
            wrk[1].resize(v.getNumSubs(), sh_order);
            wrk[2].resize(v.getNumSubs(), sh_order/2);

            Resample(   *u1, sht_upsample_, sht_, wrk[0], wrk[1], wrk[2]);
            Resample(wrk[2], sht_, sht_upsample_, wrk[0], wrk[1],    *u1);
          }
        }
        axpy(static_cast<value_type>(-1), Surf->getPosition(), *u1, *u1);
        Surf->mapToTangentSpace(*u1, false /* upsample */);

        //checks
        GeometricDot(*u1,*u1,*wrk);
        vel = MaxAbs(*wrk);
        vel = std::sqrt(vel);

        //Advecting tension (useless for implicit)
        if (params_.scheme != GloballyImplicit){
            if (params_.rep_upsample)
                WARN("Reparametrizaition is not advecting the tension in the upsample mode (fix!)");
            else {
                Surf->grad(tension_, *u2);
                GeometricDot(*u2, *u1, *wrk);
                axpy(ts/vel, *wrk, tension_, tension_);
            }
        }

        //updating position
        axpy(ts/vel, *u1, Surf->getPosition(), Surf->getPositionModifiable());

        COUTDEBUG("Iteration = "<<ii<<", |vel| = "<<vel);
    }
    INFO("Iterations = "<<ii<<", |vel| = "<<vel);

    ii=-1;
    vel=(params_.rep_tol+1);
    last_vel=(vel+1);
    for(value_type ts_=ts;ts_>1e-8;ts_*=0.5)
    while ( ++ii < params_.rep_maxit )
    {
        Surf->getSmoothedShapePosition(*u1);
        axpy(static_cast<value_type>(-1), Surf->getPosition(), *u1, *u1);

        Surf->mapToTangentSpace(*u1, false /* upsample */);

        //checks
        GeometricDot(*u1,*u1,*wrk);
        vel = MaxAbs(*wrk);
        vel = std::sqrt(vel);

        //Advecting tension (useless for implicit)
        if (params_.scheme != GloballyImplicit){
            if (params_.rep_upsample)
                WARN("Reparametrizaition is not advecting the tension in the upsample mode (fix!)");
            else {
                Surf->grad(tension_, *u2);
                GeometricDot(*u2, *u1, *wrk);
                axpy(ts_/vel, *wrk, tension_, tension_);
            }
        }

        //updating position
        axpy(ts_/vel, *u1, Surf->getPosition(), Surf->getPositionModifiable());

        COUTDEBUG("Iteration = "<<ii<<", |vel| = "<<vel);
        if (vel<params_.rep_tol || last_vel < vel) {
            flag=0;
            //if (last_vel < vel) WARN("Residual is increasing, stopping");
            break;
        }

        last_vel = vel;
    }

    { // print log(coeff)
      std::auto_ptr<Vec_t> shc = checkoutVec();
      std::auto_ptr<Vec_t> w   = checkoutVec();
      shc->replicate(Surf->getPosition());
      w  ->replicate(Surf->getPosition());
      sht_upsample_.forward(Surf->getPosition(), *w, *shc);

      {
          size_t p=shc->getShOrder();
          int ns = shc->getNumSubFuncs();
          std::vector<value_type> coeff_norm(p+1,0);
          for(int ii=0; ii<= p; ++ii){
              value_type* inPtr = shc->begin() + ii;

              int len = 2*ii + 1 - (ii/p);
              for(int jj=0; jj< len; ++jj){

                  int dist = (p + 1 - (jj + 1)/2);
                  for(int ss=0; ss<ns; ++ss){
                      coeff_norm[ii] += (*inPtr)*(*inPtr);
                      inPtr += dist;
                  }
                  inPtr--;
                  inPtr += jj%2;
              }
          }

          std::cout<<"Sph-Coeff: ";
          for(int ii=0; ii<= p; ++ii){
            std::cout<<-(int)(0.5-log(coeff_norm[ii])/log(10.0))<<' ';
          }
          std::cout<<'\n';
      }

      recycle(shc);
      recycle(w);
    }

    if (params_.rep_upsample)
        Resample(Surf->getPosition(), sht_upsample_, sht_, *u1, *u2,
            S_.getPositionModifiable());

    INFO("Iterations = "<<ii<<", |vel| = "<<vel);

    recycle(u1);
    recycle(u2);
    recycle(wrk);
    PROFILEEND("",0);

    return (flag) ? ErrorEvent::DivergenceError : ErrorEvent::Success;
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
	axpy(static_cast<value_type>(1), pos_vel_, SFb, pos_vel_);
    // } else {
    // 	axpy(static_cast<value_type>(1), SFb, pos_vel_);
    // }

    axpy(static_cast<value_type>(-1), vel, pos_vel_, vel);

    value_type err = MaxAbs(vel);
    axpy(static_cast<value_type>(1), pos_vel_, vel);

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

    Vec_t dx; dx.replicate(S_.getPosition());
    updateJacobiExplicit(S_, dt_, dx);

    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    axpy(static_cast<value_type>(-1),               dx, xnew, xnew);
    value_type err = MaxAbs(xnew);
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);
    axpy(static_cast<value_type>(1),               dx, xnew);

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
    axpy(static_cast<value_type>(1), *u2, pos_vel_, *u1);

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

    Vec_t dx; dx.replicate(S_.getPosition());
    updateJacobiGaussSeidel(S_, dt_, dx);

    axpy(static_cast<value_type>(-1), S_.getPosition(), xnew, xnew);
    axpy(static_cast<value_type>(-1),               dx, xnew, xnew);
    value_type err = MaxAbs(xnew);
    axpy(static_cast<value_type>(1), S_.getPosition(), xnew);
    axpy(static_cast<value_type>(1),               dx, xnew);

    ASSERT(err<tol, "Implicit update benchmark failed: "
        <<"err="<<err<<", tol="<<tol);
    COUT(emph<<"Implicit update benchmark passed (err="<<err<<"<tol="<<tol<<")"<<emph);

    return (err < tol);
}

#endif //NDEBUG
