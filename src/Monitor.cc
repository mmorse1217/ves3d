template<typename EvolveSurface>
MonitorBase<EvolveSurface>::~MonitorBase() 
{}

/////////////////////////////////////////////////////////////////////////////////////
template<typename EvolveSurface>
Monitor<EvolveSurface>::Monitor(const Parameters<value_type> &params) : 
    time_hor_(params.time_horizon),
    save_flag_(params.save_data),
    save_stride_(params.save_stride),
    IO(params.save_file_name),
    A0(-1),
    V0(-1)
{}

template<typename EvolveSurface>
Monitor<EvolveSurface>::~Monitor()
{}
    
template<typename EvolveSurface>
bool Monitor<EvolveSurface>::operator()(const EvolveSurface *state, 
    value_type &t, value_type &dt) 
{ 
    typename EvolveSurface::Sca_t area, vol;

    area.replicate(state->S_->getPosition());
    vol.replicate(state->S_->getPosition());
    state->S_->area(area);
    state->S_->volume(vol);
    
    value_type A(area.getDevice().MaxAbs(area.begin(), 
            state->S_->getPosition().getNumSubs()));
    value_type V( vol.getDevice().MaxAbs( vol.begin(), 
            state->S_->getPosition().getNumSubs()));
    
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
        
        if(save_flag_ && ( (static_cast<int>(t/save_stride_ - 1e-4) + 1) == static_cast<int>(t/save_stride_)))
        {
            COUT("\n           Writing data to file."<<endl);
            IO.Append(state->S_->getPosition());
            
        }
        COUT(" ------------------------------------"<<endl);
    }
    return(t<time_hor_);
}
