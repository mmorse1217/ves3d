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
    T t(0);
    T dt(params_->ts);
    T time_horizon(params_->time_horizon);

    delete F_;
    F_ = new IntVel_t(*S_, *interaction_, mats_, *params_, *vInf_, parallel_solver_);

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

    CHK( (*monitor_)( this, 0, dt) );
    INFO("Stepping with "<<params_->scheme);
    while ( ERRORSTATUS() && t < time_horizon )
    {
        CHK( (F_->*updater)(dt) );
        F_->reparam();
        t += dt;

        (*repartition_)(S_->getPositionModifiable(), F_->tension());
        CHK( (*monitor_)( this, t, dt) );
    }

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
