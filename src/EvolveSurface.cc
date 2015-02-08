template<typename T, typename DT, const DT &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t &params,
    Mats_t &mats, Vec_t &x0, BgFlow_t *vInf,  Monitor_t *M,
    Interaction_t *I, Repartition_t *R, PSolver_t *parallel_solver):
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
    S_ = new Sur_t(x0, mats_, params_.rep_up_freq, params_.rep_filter_freq);
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
}

// template<typename T, typename DT, const DT &DEVICE,
//          typename Interact, typename Repart>
// EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t &params, Mats_t &mats,
//     Sur_t *S, BgFlow_t *vInf, Monitor_t *M, Interaction_t *I, Repartition_t *R, void* user_ptr) :
//     params_(params),
//     mats_(mats),
//     S_(S),
//     vInf_(vInf),
//     monitor_(M),
//     interaction_(I),
//     repartition_(R),
//     user_ptr_(user_ptr),
//     F_(NULL),
//     ownsObjOnHeap_(false)
// {}

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
    T dt(params_.ts);
    T time_horizon(params_.time_horizon);

    delete F_;
    F_ = new IntVel_t(*S_, *interaction_, mats_, params_, *vInf_, parallel_solver_);

    //Deciding on the updater type
    Scheme_t updater(NULL);
    switch ( params_.scheme )
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
	    ErrorEvent::InvalidParameterError;
    }

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
