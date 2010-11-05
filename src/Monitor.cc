template<typename EvolveSurface>
MonitorBase<EvolveSurface>::~MonitorBase() 
{}

/////////////////////////////////////////////////////////////////////////////////////
template<typename EvolveSurface>
Monitor<EvolveSurface>::Monitor(const Parameters<value_type> &params) : 
    save_flag_(params.save_data),
    save_stride_(params.save_stride),
    IO(params.save_file_name, 5 * (DIM + 1) * /* position + tension */ 
        gridDimOf(params.sh_order).first * 
        gridDimOf(params.sh_order).second * sizeof(value_type)),
    A0(-1),
    V0(-1),
    last_save(-1),
    params_(params)
{}

template<typename EvolveSurface>
Monitor<EvolveSurface>::~Monitor()
{
    //Due to a bug'ish in the DataIO class
    IO.FlushBuffer<value_type>();
}
    
template<typename EvolveSurface>
Error_t Monitor<EvolveSurface>::operator()(const EvolveSurface *state, 
    const value_type &t, value_type &dt) 
{
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
        
        bool save_now = static_cast<int>(t/save_stride_) > last_save;

        if ( save_flag_ && save_now )
        {
            COUT("\n           Writing data to file."<<endl);
            IO.Append(state->S_->getPosition());
            IO.Append(state->F_->tension());
            last_save++;
        }
        COUT(" ------------------------------------"<<endl);
    }
    
    Error_t return_val(Success);
    if ( abs(A/A0-1) > params_.error_factor  || abs(V/V0-1) > params_.error_factor )
        return_val = AccuracyError;

    return return_val;
}
