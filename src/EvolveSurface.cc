template<typename SurfContainer>
Monitor<SurfContainer>::Monitor(const Parameters<value_type> &params) : 
    time_hor_(params.time_horizon),
    save_flag_(params.save_data),
    save_stride_(params.save_stride),
    IO(params.save_file_name),
    A0(-1),
    V0(-1)
{};
    
template<typename SurfContainer>
bool Monitor<SurfContainer>::operator()(const SurfContainer &state, 
    value_type &t, value_type &dt, 
    void*/* to be compatible with the general signature, not used*/)
{
    ///@todo move temporary object
    typename SurfContainer::Sca_t area, vol;
    area.replicate(state.getPosition());
    vol.replicate(state.getPosition());
    state.area(area);
    state.volume(vol);
    
    value_type A(area.getDevice().MaxAbs(area.begin(), 
            state.getPosition().getNumSubs()));
    value_type V( vol.getDevice().MaxAbs( vol.begin(), 
            state.getPosition().getNumSubs()));
    
    if(A0 == -1)
    {
        A0 = A;
        V0 = V;
    }
    
#pragma omp critical (monitor)
    {
        COUT("\n  Monitor :"
            <<"\n           thread       = "<<omp_get_thread_num()
            <<"/"<<omp_get_num_threads()
            <<"\n           t            = "<<fixed<<t
            <<scientific<<setprecision(4)
            <<"\n           area   error = "<<abs(A/A0-1)
            <<scientific<<setprecision(4)
            <<"\n           volume error = "<<abs(V/V0-1)<<endl);
        
        if(save_flag_ && (t/save_stride_ == static_cast<int>(t/save_stride_)))
        {
            COUT("\n           Writing data to file."<<endl);
            IO.Append(state.getPosition());
            
        }
        COUT(" ------------------------------------"<<endl);
    }
    return(t<time_hor_);
}

template<typename Container, 
         typename Interaction,
         typename Mntr,
         typename Repart>
EvolveSurface<Container, Interaction, Mntr, Repart>::EvolveSurface(
    OperatorsMats<typename Container::Sca_t> &mats, 
    const Parameters<value_type> &params, Mntr& monitor, 
    Repart& repartion) : 
    mats_(mats), 
    params_(params), 
    monitor_(monitor), 
    repartition_(repartion)
{}

template<typename Container, 
         typename Interaction,
         typename Mntr,
         typename Repart>
void EvolveSurface<Container, Interaction, Mntr, Repart>::operator()(
    Container &S, Interaction &Inter, void* usr_ptr)
{
    value_type t(0);
    value_type dt(params_.ts);
    IntVel_t F(S, Inter, mats_, params_);
    F.usr_ptr_ = usr_ptr;
        
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
    
    while ( monitor_(S, t, dt, usr_ptr) )
    {
        (F.*updater)(dt);       
        F.reparam();
        t += dt;
        
        // repartition_.operator()<Container::Sca_t>(S.getPositionModifiable(), 
        //     F.tension(), usr_ptr);
    }
}
