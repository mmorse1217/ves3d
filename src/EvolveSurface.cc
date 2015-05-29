template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t *params,
    const Mats_t &mats, BgFlow_t *vInf,  Monitor_t *M,
    Interaction_t *I, Repartition_t *R, PSolver_t *parallel_solver, Vec_t *x0):
    params_(params),
    mats_(mats),
    vInf_(vInf),
    parallel_solver_(parallel_solver),
    monitor_(M),
    interaction_(I),
    repartition_(R),
    F_(NULL)
{
    objsOnHeap_[0] = objsOnHeap_[1] = objsOnHeap_[2] = false;

    S_ = new Sur_t(params->sh_order, mats_, x0, params_->filter_freq, params_->rep_filter_freq);
    S_->set_name("surface");

    if ( monitor_ == NULL)
    {
        monitor_ = new Monitor<EvolveSurface>(params_);
        objsOnHeap_[0] = true;
    }

    if ( interaction_ == NULL)
    {
        interaction_ = new Interaction_t();
        objsOnHeap_[1] = true;
    }

    if ( repartition_ == NULL)
    {
        repartition_ = new Repartition_t();
        objsOnHeap_[2] = true;
    }
    INFO("Created a new object");
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::~EvolveSurface()
{
    delete S_;
    delete F_;

    if ( objsOnHeap_[0] ) delete monitor_;
    if ( objsOnHeap_[1] ) delete interaction_;
    if ( objsOnHeap_[2] ) delete repartition_;
}

template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::Evolve()
{
    PROFILESTART();
    T t(0);
    T dt(params_->ts);
    T time_horizon(params_->time_horizon);

    delete F_;
    F_ = new IntVel_t(*S_, *interaction_, mats_, *params_, *vInf_, parallel_solver_);
    IntVel_t* F_coarse_ = new IntVel_t(*S_, *interaction_, mats_, *params_, *vInf_, parallel_solver_);

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
	    return ErrorEvent::NotImplementedError;
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
    TimeAdaptive time_adap=(params_->time_adaptive?TimeAdapErrAreaVol:TimeAdapNone);

    Vec_t dx, x0, x_coarse;
    CHK( (*monitor_)( this, 0, dt) );
    INFO("Stepping with "<<params_->scheme);
    while ( ERRORSTATUS() && t < time_horizon )
    {
        if(time_adap==TimeAdapErr){ // Adaptive using 2*dt time-step for error
            Error_t err=ErrorEvent::Success;

            // Copy S_->getPosition
            x0.replicate(S_->getPosition());
            axpy(static_cast<value_type>(0.0), S_->getPosition(), S_->getPosition(), x0);

            // 2*dt time-step
            x_coarse.replicate(S_->getPosition());
            if(err==ErrorEvent::Success) err=(F_coarse_->*updater)(*S_, 2*dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), x_coarse);

            // 2 x (dt time-step)
            if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());
            if(err==ErrorEvent::Success) err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());

            // Compute error
            axpy(static_cast<value_type>(-1.0), S_->getPosition(), x_coarse, x_coarse);
            //value_type error=MaxAbs(x_coarse);
            value_type error=sqrt(AlgebraicDot(x_coarse,x_coarse)/x_coarse.size());

            bool accept=true;
            value_type dt_new=dt;
            { // Compute dt_new
                value_type timestep_order=1;
                value_type time_horizon=params_->time_horizon;

                value_type beta;
                beta = (1.0/error) * (dt/time_horizon) * params_->error_factor;
                beta = std::pow<value_type>(beta, timestep_order);
                if(err!=ErrorEvent::Success) beta=0.5;

                beta=std::min(beta,1.5);
                beta=std::max(beta,0.5);

                value_type beta_scale=std::pow(sqrt(0.9),timestep_order) * beta;
                if(beta_scale<0.5) accept=false;
                dt_new=beta_scale * dt;
            }
            if(accept){ // Increment t
                t += 2*dt;
            }else{ // Restore original S_
                axpy(static_cast<value_type>(0.0), x0, x0, S_->getPositionModifiable());
            }

            INFO("Time-adaptive: error/dt = "<<error/dt<<", dt_new = "<<dt_new);
            dt=dt_new;
        }else if(time_adap==TimeAdapErrAreaVol){ // Adaptive using area, volume error
            Error_t err=ErrorEvent::Success;

            // Copy S_->getPosition
            x0.replicate(S_->getPosition());
            axpy(static_cast<value_type>(0.0), S_->getPosition(), S_->getPosition(), x0);

            // Compute initial area/volume
            value_type A0, V0;
            static Sca_t area0, vol0;
            size_t N_ves=S_->getPosition().getNumSubs();
            area0.replicate(S_->getPosition()); S_->area  (area0);
            vol0 .replicate(S_->getPosition()); S_->volume( vol0);
            A0=Sca_t::getDevice().MaxAbs(area0.begin(), N_ves);
            V0=Sca_t::getDevice().MaxAbs( vol0.begin(), N_ves);

            // dt time-step
            err=(F_->*updater)(*S_, dt, dx);
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());

            // Compute area/volume error
            value_type A_err, V_err;
            static Sca_t area_err, vol_err;
            area_err.replicate(S_->getPosition()); S_->area  (area_err);
            vol_err .replicate(S_->getPosition()); S_->volume( vol_err);
            axpy(static_cast<value_type>(-1.0), area0, area_err, area_err);
            axpy(static_cast<value_type>(-1.0),  vol0,  vol_err,  vol_err);
            A_err=Sca_t::getDevice().MaxAbs(area_err.begin(), N_ves);
            V_err=Sca_t::getDevice().MaxAbs( vol_err.begin(), N_ves);

            bool accept=true;
            value_type dt_new=dt;
            { // Compute dt_new
                value_type timestep_order=1;
                value_type time_horizon=params_->time_horizon;

                value_type beta_A;
                beta_A = (A0/A_err) * (dt/time_horizon) * params_->error_factor;
                beta_A = std::pow<value_type>(beta_A, timestep_order);
                if(err!=ErrorEvent::Success) beta_A=0.5;

                value_type beta_V;
                beta_V = (V0/V_err) * (dt/time_horizon) * params_->error_factor;
                beta_V = std::pow<value_type>(beta_V, timestep_order);
                if(err!=ErrorEvent::Success) beta_V=0.5;

                value_type beta=std::min(beta_A,beta_V);
                beta=std::min(beta,1.5);
                beta=std::max(beta,0.5);

                value_type beta_scale=std::pow(sqrt(0.9),timestep_order) * beta;
                if(beta_scale<0.5) accept=false;
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
            CHK( (F_->*updater)(*S_, dt, dx) );
            axpy(static_cast<value_type>(1.0), dx, S_->getPosition(), S_->getPositionModifiable());

            t += dt;
        }

        F_->reparam();
        (*repartition_)(S_->getPositionModifiable(), F_->tension());
        CHK( (*monitor_)( this, t, dt) );
    }

    delete F_coarse_;

    PROFILEEND("",0);
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
    int ii;

    is>>s;
    ASSERT(s=="EVOLVE", "Bad input string (missing header).");

    is>>key>>s;
    ASSERT(key=="version:", "bad key version");
    INFO("Unpacking checkpoint from version "<<s<<" (current version "<<VERSION<<")");

    is>>key>>Streamable::name_;
    ASSERT(key=="name:", "bad key name");
    S_->unpack(is,format);
    is>>s;
    ASSERT(s=="/EVOLVE", "Bad input string (missing footer).");

    return ErrorEvent::Success;
}
