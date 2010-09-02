template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t &params, Mats_t &mats, 
    Vec_t &x0, BgFlow_t *vInf, InteractionFun_t interaction_handle, void* user_ptr,
    GlobalRepart_t repartition_handle) :
    params_(params),
    mats_(mats),
    vInf_(vInf),
    user_ptr_(user_ptr),
    ownsObjOnHeap_(false)
{
    ownsObjOnHeap_ = true;
    S_ = new Sur_t(x0, mats_);
    monitor_ = new Monitor<EvolveSurface>(params_);
    interaction_ = new Interaction_t(interaction_handle);
    repartition_ = new Repartition_t(repartition_handle);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::EvolveSurface(Params_t &params, Mats_t &mats, 
    Sur_t *S, BgFlow_t *vInf, Monitor_t *M, Interaction_t *I, Repartition_t *R, void* user_ptr) :
    params_(params),
    mats_(mats),
    S_(S),
    vInf_(vInf),
    monitor_(M),
    interaction_(I),
    repartition_(R),
    user_ptr_(user_ptr),
    ownsObjOnHeap_(false)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
EvolveSurface<T, DT, DEVICE, Interact, Repart>::~EvolveSurface()
{
    if ( ownsObjOnHeap_ )
    {
        delete S_;
        delete monitor_;
        delete interaction_;
        delete repartition_;
    }
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE,
         typename Interact, typename Repart>
MonitorReturn EvolveSurface<T, DT, DEVICE, Interact, Repart>::Evolve() 
{
    T t(0);
    T dt(params_.ts);
    IntVel_t F(*S_, *interaction_, mats_, params_, *vInf_);
    F.usr_ptr_ = user_ptr_;
     
    //Deciding on the updater type
    Scheme_t updater;
    switch ( params_.scheme )
    {
        case Explicit:
            updater = &IntVel_t::updatePositionExplicit; 
            break;

        case SemiImplicit:
            updater = &IntVel_t::updatePositionImplicit; 
            break;
    }
    
    enum MonitorReturn monitor_status(StatusOK);

    while ( monitor_status == StatusOK)
    {
        (F.*updater)(dt);       
        F.reparam();
        t += dt;
        
        //repartition_.operator()<Container::Sca_t>(S.getPositionModifiable(), 
        //    F.tension(), usr_ptr);
        
        monitor_status = (*monitor_)( this, t, dt);
    }
    return monitor_status;
}
