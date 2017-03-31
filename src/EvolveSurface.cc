template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t *params,
    const Mats_t &mats, BgFlow_t *vInf,  Monitor_t *M, Interaction_t *I,
    Repartition_t *R, PSolver_t *parallel_solver, Vec_t *x0,VProp_t *ves_props):
    params_(params),
    ves_props_(ves_props),
    mats_(mats),
    S_up_(NULL),
    vInf_(vInf),
    parallel_solver_(parallel_solver),
    monitor_(M),
    interaction_(I),
    repartition_(R),
    F_(NULL)
{
    ownedObjs_[0] = ownedObjs_[1] = ownedObjs_[2] = ownedObjs_[3] = false;

    S_ = new Sur_t(params->sh_order, mats_, x0, params_->filter_freq,
        params_->rep_filter_freq,params_->rep_type,params_->rep_exponent);
    S_->set_name("surface");

    if ( monitor_ == NULL)
    {
        monitor_ = new Monitor<EvolveSurface>(params_);
        ownedObjs_[0] = true;
    }

    if ( interaction_ == NULL)
    {
        interaction_ = new Interaction_t();
        ownedObjs_[1] = true;
    }

    if ( repartition_ == NULL)
    {
        repartition_ = new Repartition_t();
        ownedObjs_[2] = true;
    }

    if (ves_props_ == NULL)
    {
        ves_props_ = new VProp_t();
        ves_props_->setFromParams(*params_);
        ownedObjs_[3] = true;
    }

    set_name("evolve_surface");
    INFO("Created a new object");
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::~EvolveSurface()
{
    delete S_;
    delete F_;

    if ( ownedObjs_[0] ) delete monitor_;
    if ( ownedObjs_[1] ) delete interaction_;
    if ( ownedObjs_[2] ) delete repartition_;
    if ( ownedObjs_[3] ) delete ves_props_;
    if (S_up_) delete S_up_;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::ReinitInterfacialVelocity()
{
    delete F_;
    F_ = new IntVel_t(*S_, *interaction_, mats_, *params_, *ves_props_,
        *vInf_, parallel_solver_);
    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::Evolve()
{
    PROFILESTART();
    T t(0);
    T dt(params_->ts);
    T time_horizon(params_->time_horizon);

    INFO("The vesicles' material properties:\n"<<ves_props_);

    ReinitInterfacialVelocity();
    //Deciding on the updater type
    Scheme_t updater(NULL);
    switch ( params_->scheme )
    {
	case JacobiBlockExplicit:
            updater = &IntVel_t::updateJacobiExplicit;
            break;

	case JacobiBlockGaussSeidel:
            updater = &IntVel_t::updateJacobiGaussSeidel;
            break;

	case JacobiBlockImplicit:
	    updater = &IntVel_t::updateJacobiImplicit;
	    break;

	case GloballyImplicit:
	    updater = &IntVel_t::updateImplicit;
	    break;

	default:
	  return ErrorEvent::InvalidParameterError;
    }

    enum TimeAdaptive{
      TimeAdapErr,
      TimeAdapErrAreaVol,
      TimeAdapNone
    };
    TimeAdaptive time_adap=(params_->time_adaptive?TimeAdapErr:TimeAdapNone);

    Sca_t area, vol;
    { // Compute area, vol
        S_->resample(params_->upsample_freq, &S_up_); // up-sample
        int N_ves=S_up_->getNumberOfSurfaces();
        area.resize(N_ves,1); S_up_->area  (area);
        vol .resize(N_ves,1); S_up_->volume( vol);
        //@bug downsample seems unnecessary
        //S_up_->resample(params_->sh_order, &S_); // down-sample
    }

    { // Sanity check for time-stepper
        Vec_t y0, y1;
        Error_t err=ErrorEvent::Success;
        if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, y0);
        value_type max_y0=MaxAbs(y0);
        if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, y1);
        axpy(static_cast<value_type>(-1.0), y0, y1, y0);
        value_type max_err=MaxAbs(y0);

        { // max_err = MPI_MAX(max_err)
            ASSERT(typeid(T)==typeid(double),"Only works for double"); // @bug this only works for T==double
            value_type loc;

            loc=max_y0;
            MPI_Allreduce(&loc, &max_y0, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);
            loc=max_err;
            MPI_Allreduce(&loc, &max_err, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);
        }
        //@bug I (ABT) think tolerance is true when solving for
        //velocity, not for position
        INFO("Sanity check error="<<max_err<<
            ", relative error="<<max_err/max_y0<<
            ", rtol="<<dt*params_->time_tol);

        //if(max_err>max_y0*dt*params_->time_tol)
        if(max_err>max_y0*dt*params_->time_tol && params_->scheme==GloballyImplicit)
            CERR_LOC("Sanity check for time-stepper failed!", std::endl,exit(1));
    }

    Vec_t dx, x0, x_dt, x_2dt;
    CHK( (*monitor_)( this, 0, dt) );
    INFO("Stepping with "<<params_->scheme);

    MPI_Comm comm=MPI_COMM_WORLD;
    pvfmm::Profile::Enable(true);
    while ( ERRORSTATUS() && t < time_horizon && dt>1e-10 )
    {
        pvfmm::Profile::Tic("TimeStep",&comm,true);

        if(time_adap==TimeAdapErr){ // Adaptive using 2*dt time-step for error
            dt=std::min((time_horizon-t)/2, dt);
            Error_t err=ErrorEvent::Success;

            // Copy S_->getPosition
            x0.replicate(S_->getPosition());
            axpy(static_cast<value_type>(0.0), S_->getPosition(), S_->getPosition(), x0);

            // dt time-step
            pvfmm::Profile::Tic("GMRES1",&comm,true);
            x_dt.replicate(S_->getPosition());
            if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), x_dt);
            pvfmm::Profile::Toc();

            // Check integration error
            value_type stokes_error=F_->StokesError(x_dt);

            // 2*dt time-step
            pvfmm::Profile::Tic("GMRES2",&comm,true);
            x_2dt.replicate(S_->getPosition());
            axpy(static_cast<value_type>(0.0), x0, x0, S_->getPositionModifiable());
            if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, 2*dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), x_2dt);
            pvfmm::Profile::Toc();

            // dt time-step
            pvfmm::Profile::Tic("GMRES3",&comm,true);
            axpy(static_cast<value_type>(0.0), x_dt, x_dt, S_->getPositionModifiable());
            if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());
            pvfmm::Profile::Toc();

            // Check integration error
            stokes_error=std::max(stokes_error, F_->StokesError(S_->getPosition()));
            { // stokes_error = MPI_MAX(stokes_error)
              assert(typeid(T)==typeid(double)); // @bug this only works for T==double
              value_type error_loc=stokes_error;
              MPI_Allreduce(&error_loc, &stokes_error, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);
            }

            // Compute error
            axpy(static_cast<value_type>(-1.0), S_->getPosition(), x_2dt, x_2dt);
            value_type error=MaxAbs(x_2dt);
            { // error = MPI_MAX(error)
              assert(typeid(T)==typeid(double)); // @bug this only works for T==double
              value_type error_loc=error;
              MPI_Allreduce(&error_loc, &error, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);
            }

            int accept=1;
            value_type dt_new=dt;
            { // Compute dt_new
                value_type timestep_order=1;
                value_type time_horizon=params_->time_horizon;

                value_type beta;
                beta = (1.0/error) * (dt/time_horizon) * params_->error_factor;
                beta = std::pow(beta, timestep_order);
                if(err!=ErrorEvent::Success) beta=0.5;
                if(stokes_error*dt>params_->time_tol) beta=params_->time_tol/(stokes_error*dt); // This is required for GMRES to converge

                beta=std::min(beta,1.5);
                beta=std::max(beta,0.5);

                value_type beta_scale=std::pow(sqrt(0.9),timestep_order) * beta;
                accept=(beta_scale<0.5?0:1);
                dt_new=beta_scale * dt;
            }
            if(accept){ // Increment t
                t += 2*dt;
            }else{ // Restore original S_
                axpy(static_cast<value_type>(0.0), x0, x0, S_->getPositionModifiable());
            }

            INFO("Time-adaptive: error/dt = "<<error/dt<<", error/dt^2 = "<<error/dt/dt<<", dt_new = "<<dt_new);
            dt=dt_new;
        }else if(time_adap==TimeAdapErrAreaVol){ // Adaptive using area, volume error
            dt=std::min(time_horizon-t, dt);
            Error_t err=ErrorEvent::Success;

            // Copy S_->getPosition
            x0.replicate(S_->getPosition());
            axpy(static_cast<value_type>(0.0), S_->getPosition(), S_->getPosition(), x0);

            // Compute initial area/volume
            value_type A0, V0;
            static Sca_t area0, vol0;
            size_t N_ves=S_->getPosition().getNumSubs();
            S_->resample(params_->upsample_freq, &S_up_); // up-sample
            area0.replicate(S_up_->getPosition()); S_up_->area  (area0);
            vol0 .replicate(S_up_->getPosition()); S_up_->volume( vol0);
            A0=Sca_t::getDevice().MaxAbs(area0.begin(), N_ves);
            V0=Sca_t::getDevice().MaxAbs( vol0.begin(), N_ves);

            // dt time-step
            pvfmm::Profile::Tic("GMRES",&comm,true);
            err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());
            pvfmm::Profile::Toc();

            // Compute area/volume error
            value_type A_err, V_err;
            static Sca_t area_err, vol_err;
            S_->resample(params_->upsample_freq, &S_up_); // up-sample
            area_err.replicate(S_up_->getPosition()); S_up_->area  (area_err);
            vol_err .replicate(S_up_->getPosition()); S_up_->volume( vol_err);
            axpy(static_cast<value_type>(-1.0), area0, area_err, area_err);
            axpy(static_cast<value_type>(-1.0),  vol0,  vol_err,  vol_err);
            A_err=Sca_t::getDevice().MaxAbs(area_err.begin(), N_ves);
            V_err=Sca_t::getDevice().MaxAbs( vol_err.begin(), N_ves);
            { // error = MPI_MAX(error)
              assert(typeid(T)==typeid(double)); // @bug this only works for T==double
              value_type error_loc;

              error_loc=A_err;
              MPI_Allreduce(&error_loc, &A_err, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);

              error_loc=V_err;
              MPI_Allreduce(&error_loc, &V_err, 1, MPI_DOUBLE, MPI_MAX, VES3D_COMM_WORLD);
            }

            int accept=1;
            value_type dt_new=dt;
            { // Compute dt_new
                value_type timestep_order=1;
                value_type time_horizon=params_->time_horizon;

                value_type beta_A;
                beta_A = (A0/A_err) * (dt/time_horizon) * params_->error_factor;
                beta_A = std::pow(beta_A, timestep_order);
                if(err!=ErrorEvent::Success) beta_A=0.5;

                value_type beta_V;
                beta_V = (V0/V_err) * (dt/time_horizon) * params_->error_factor;
                beta_V = std::pow(beta_V, timestep_order);
                if(err!=ErrorEvent::Success) beta_V=0.5;

                value_type beta=std::min(beta_A,beta_V);
                beta=std::min(beta,1.5);
                beta=std::max(beta,0.5);

                value_type beta_scale=std::pow(sqrt(0.9),timestep_order) * beta;
                accept=(beta_scale<0.5?0:1);
                dt_new=beta_scale * dt;
            }
            if(accept){ // Increment t
                t += dt;
            }else{ // Restore original S_
                axpy(static_cast<value_type>(0.0), x0, x0, S_->getPositionModifiable());
            }

            INFO("Time-adaptive: A_err/dt = "<<(A_err/A0)/dt<<", V_err/dt = "<<(V_err/V0)/dt<<", dt_new = "<<dt_new);
            dt=dt_new;
        }else if(time_adap==TimeAdapNone){ // No adaptive
            pvfmm::Profile::Tic("GMRES",&comm,true);
            CHK( (F_->*updater)(*S_, dt, dx) );
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());
            pvfmm::Profile::Toc();

            t += dt;
        }

        pvfmm::Profile::Tic("Reparam",&comm,true);
        F_->reparam();
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("AreaVolume",&comm,true);
        AreaVolumeCorrection(area, vol);
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("Repartition",&comm,true);
        (*repartition_)(S_->getPositionModifiable(), F_->tension());
        pvfmm::Profile::Toc();
        pvfmm::Profile::Tic("Monitor",&comm,true);
        CHK( (*monitor_)( this, t, dt) );
        pvfmm::Profile::Toc();

        pvfmm::Profile::Toc();
        //pvfmm::Profile::print(&comm);
    }
    PROFILEEND("",0);
    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::AreaVolumeCorrection(const Sca_t& area, const Sca_t& vol, const value_type tol)
{
    static GaussLegendreIntegrator<Sca_t> integrator;
    const DT& device=Sca_t::getDevice();
    int N_ves=S_->getNumberOfSurfaces();
    int iter(-1);

    // store collision free state in x_old
    Vec_t x_old, u1;
    x_old.replicate(S_->getPosition());
    u1.replicate(S_->getPosition());
    axpy(static_cast<value_type>(1.0), S_->getPosition(), x_old);
    // end of store 
    
    while (++iter < params_->rep_maxit){
        S_->resample(params_->upsample_freq, &S_up_); // up-sample

        const Vec_t& Normal  =S_up_->getNormal();
        const Sca_t& AreaElem=S_up_->getAreaElement();
        const Sca_t& MeanCurv=S_up_->getMeanCurv();

        Sca_t X; // First perturbation direction
        { // X = -2.0*MeanCurv
            X.replicate(MeanCurv);
            axpy(-2.0, MeanCurv, X);
        }

        Sca_t Y; // Second perturbation direction
        { // Y = [1, ..., 1]
            Y.replicate(MeanCurv);
            device.Memset(Y.begin(), 1, Y.size()*sizeof(value_type));
            xyInv(Y,Y,Y);
        }

        Sca_t dX(N_ves,1);
        Sca_t dY(N_ves,1);
        { // Set dX, dY
            Sca_t area_err, vol_err;
            { // compute error
                area_err.resize(N_ves,1); S_up_->area  (area_err);
                vol_err .resize(N_ves,1); S_up_->volume( vol_err);
                device.axpy(-1.0, area_err.begin(), area.begin(), N_ves, area_err.begin());
                device.axpy(-1.0,  vol_err.begin(),  vol.begin(), N_ves,  vol_err.begin());

                value_type area_max_err=device.MaxAbs<value_type>(area_err.begin(),N_ves)/device.MaxAbs<value_type>(area.begin(),N_ves);
                value_type  vol_max_err=device.MaxAbs<value_type>( vol_err.begin(),N_ves)/device.MaxAbs<value_type>( vol.begin(),N_ves);
                COUTDEBUG("Iteration = "<<iter<<", area error="<<area_max_err<<", vol error="<<vol_max_err);
                if(std::max(area_max_err, vol_max_err)<tol) break;
            }

            Sca_t dAdX(N_ves,1);
            Sca_t dAdY(N_ves,1);
            Sca_t dVdX(N_ves,1);
            Sca_t dVdY(N_ves,1);
            { // dA/dX
                Sca_t tmp; tmp.replicate(MeanCurv);
                xy(X, X, tmp); integrator(tmp,AreaElem,dAdX);
            }
            { // dA/dY
                Sca_t tmp; tmp.replicate(MeanCurv);
                xy(X, Y, tmp); integrator(tmp,AreaElem,dAdY);
            }
            { // dV/dX
                Sca_t tmp; tmp.replicate(MeanCurv);
                xy(Y, X, tmp); integrator(tmp,AreaElem,dVdX);
            }
            { // dV/dY
                Sca_t tmp; tmp.replicate(MeanCurv);
                xy(Y, Y, tmp); integrator(tmp,AreaElem,dVdY);
            }

            Sca_t DetInv(N_ves,1);
            { // Set DetInv = (dA/dX.dV/dY - dA/dY.dV/dX)^-1
                device.xy(dAdX.begin(), dVdY.begin(), N_ves, dX.begin());
                device.xy(dAdY.begin(), dVdX.begin(), N_ves, dY.begin());
                device.axpy(-1.0, dY.begin(), dX.begin(), N_ves, DetInv.begin());
                device.xyInv<value_type>(NULL, DetInv.begin(), N_ves, DetInv.begin());
            }

            Sca_t dXdA(N_ves,1);
            Sca_t dXdV(N_ves,1);
            Sca_t dYdA(N_ves,1);
            Sca_t dYdV(N_ves,1);
            { // dX/dA =  dVdY * DetInv
                device.xy(dVdY.begin(), DetInv.begin(), N_ves, dXdA.begin());
            }
            { // dX/dV = -dAdY * DetInv
                device.xy(dAdY.begin(), DetInv.begin(), N_ves, dXdV.begin());
                device.axpy<value_type>(-1.0, dXdV.begin(), NULL, N_ves, dXdV.begin());
            }
            { // dY/dA = -dVdX * DetInv
                device.xy(dVdX.begin(), DetInv.begin(), N_ves, dYdA.begin());
                device.axpy<value_type>(-1.0, dYdA.begin(), NULL, N_ves, dYdA.begin());
            }
            { // dY/dV =  dAdX * DetInv
                device.xy(dAdX.begin(), DetInv.begin(), N_ves, dYdV.begin());
            }

            { // dX = dX/dA*area_err + dX/dV*vol_err
              device.xy(dXdA.begin(), area_err.begin(), N_ves, dXdA.begin());
              device.xy(dXdV.begin(),  vol_err.begin(), N_ves, dXdV.begin());
              device.axpy(1.0, dXdA.begin(), dXdV.begin(), N_ves, dX.begin());
            }
            { // dY = dY/dA*area_err + dY/dV*vol_err
              device.xy(dYdA.begin(), area_err.begin(), N_ves, dYdA.begin());
              device.xy(dYdV.begin(),  vol_err.begin(), N_ves, dYdV.begin());
              device.axpy(1.0, dYdA.begin(), dYdV.begin(), N_ves, dY.begin());
            }
        }

        Vec_t dS; dS.replicate(Normal);
        Vec_t& position=S_up_->getPositionModifiable();
        { // position += dX*X.*Normal
            device.xvpw<value_type>( X.begin(), Normal.begin(), NULL, Normal.getStride(), Normal.getNumSubs(), dS.begin());
            device.avpw<value_type>(dX.begin(),     dS.begin(), NULL, Normal.getStride(), Normal.getNumSubs(), dS.begin());
            axpy(1, dS, position, position);
        }
        { // position += dY*Y.*Normal
            device.xvpw<value_type>( Y.begin(), Normal.begin(), NULL, Normal.getStride(), Normal.getNumSubs(), dS.begin());
            device.avpw<value_type>(dY.begin(),     dS.begin(), NULL, Normal.getStride(), Normal.getNumSubs(), dS.begin());
            axpy(1, dS, position, position);
        }

        S_up_->resample(params_->sh_order, &S_); // down-sample
    }
    INFO("Number of iterations : "<<iter);
    
    // begin for collision
    // project u1 to collision free
    INFO("Begin Project AreaVolume correction direction to without contact.");
    axpy(static_cast<value_type>(-1.0), x_old, S_->getPosition(), u1);
    F_->projectU1(u1, x_old);
    axpy(static_cast<value_type>(1.0), u1, x_old, S_->getPositionModifiable());
    INFO("End Project AreaVolume correction  direction to without contact.");
    // end for collision

    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::pack(
    std::ostream &os, Streamable::Format format) const{

    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");

    os<<"EVOLVE\n";
    os<<"version: "<<VERSION<<"\n";
    os<<"name: "<<Streamable::name_<<"\n";
    ves_props_->pack(os,format);
    S_->pack(os,format);
    os<<"/EVOLVE\n";

    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::unpack(
    std::istream &is, Streamable::Format format){

    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s,key;
    int version;

    is>>s;
    ASSERT(s=="EVOLVE", "Bad input string (missing header).");

    is>>key>>version;
    ASSERT(key=="version:", "bad key version");

    is>>key;
    if (key=="+") {++version;is>>key;}
    is>>Streamable::name_;
    ASSERT(key=="name:", "bad key name");

    if (version>590){
        ves_props_->unpack(is,format);
    } else {
        ves_props_->setFromParams(*params_);
    }
    S_->unpack(is,format);
    is>>s;
    ASSERT(s=="/EVOLVE", "Bad input string (missing footer).");

    INFO("Unpacked "<<Streamable::name_<<" data from version "<<version<<" (current version "<<VERSION<<")");
    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::getSurfaceUp(const Sur_t *&S_up) const
{
    S_->resample(params_->upsample_freq, &S_up_);
    S_up = S_up_;
    return ErrorEvent::Success;
}
