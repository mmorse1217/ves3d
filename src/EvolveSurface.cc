template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t &params,
    Mats_t &mats, Vec_t &x0, BgFlow_t *vInf,  Monitor_t *M,
    Interaction_t *I, Repartition_t *R, void* user_ptr):
    params_(params),
    mats_(mats),
    vInf_(vInf),
    monitor_(M),
    interaction_(I),
    repartition_(R),
    user_ptr_(user_ptr),
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

// template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
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

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::~EvolveSurface()
{
    delete S_;
    delete F_;

    if ( objsOnHeap_[0] ) delete monitor_;
    if ( objsOnHeap_[1] ) delete interaction_;
    if ( objsOnHeap_[2] ) delete repartition_;

}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
Error_t EvolveSurface<T, DT, DEVICE, Interact, Repart>::Evolve()
{
    T t(0);
    T dt(params_.ts);
    T time_horizon(params_.time_horizon);

    delete F_;
    F_ = new IntVel_t(*S_, *interaction_, mats_, params_, *vInf_);
    F_->usr_ptr_ = user_ptr_;

    //Deciding on the updater type
    Scheme_t updater;
    switch ( params_.scheme )
    {
        case Explicit:
            updater = &IntVel_t::updatePositionExplicit;
            break;

        case BlockImplicit:
            updater = &IntVel_t::updatePositionImplicit;
            break;
    }

    while ( ERRORSTATUS() && t < time_horizon )
    {
        QC( (F_->*updater)(dt) );
        F_->reparam();
        t += dt;

        (*repartition_)(S_->getPositionModifiable(), F_->tension(), user_ptr_);
        QC( (*monitor_)( this, t, dt) );
    }

    return Success;
}
