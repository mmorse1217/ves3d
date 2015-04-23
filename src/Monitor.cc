template<typename EvolveSurface>
MonitorBase<EvolveSurface>::~MonitorBase()
{}

/////////////////////////////////////////////////////////////////////////////////////
///@todo move the buffer size to the parameters
template<typename EvolveSurface>
Monitor<EvolveSurface>::Monitor(const Parameters<value_type> *params) :
    checkpoint_flag_(params->checkpoint),
    checkpoint_stride_(params->checkpoint_stride),
    A0_(-1),
    V0_(-1),
    last_checkpoint_(-1),
    time_idx_(-1),
    params_(params)
{}

template<typename EvolveSurface>
Monitor<EvolveSurface>::~Monitor()
{}

template<typename EvolveSurface>
Error_t Monitor<EvolveSurface>::operator()(const EvolveSurface *state,
    const value_type &t, value_type &dt)
{
    area_.replicate(state->S_->getPosition());
    vol_.replicate(state->S_->getPosition());
    state->S_->area(area_);
    state->S_->volume(vol_);

    value_type A(area_.getDevice().MaxAbs(area_.begin(),
            state->S_->getPosition().getNumSubs()));
    value_type V( vol_.getDevice().MaxAbs(vol_.begin(),
            state->S_->getPosition().getNumSubs()));

    if(A0_ == -1)
    {
        A0_ = A;
        V0_ = V;
    }

#pragma omp critical (monitor)
    {

        INFO("Monitor: thread = "<<omp_get_thread_num()<<"/"<<omp_get_num_threads()
            <<", t = "<<SCI_PRINT_FRMT<<t
            <<", area error = "<<SCI_PRINT_FRMT<<fabs(A/A0_-1)
            <<", volume error = "<<SCI_PRINT_FRMT<<fabs(V/V0_-1));

        bool checkpoint_now = static_cast<int>(t/checkpoint_stride_) > last_checkpoint_;

        if ( checkpoint_flag_ && checkpoint_now )
        {
	    ++time_idx_;
	    std::string fname(params_->checkpoint_file_name);
	    char suffix[6];
	    sprintf(suffix, "%05d", time_idx_);
	    d_["time_idx"] = std::string(suffix);
	    expand_template(&fname, d_);

	    std::stringstream ss;
	    ss<<std::scientific<<std::setprecision(16);
	    params_->pack(ss, Streamable::ASCII);
	    state->pack(ss, Streamable::ASCII);

	    INFO("Writing data to file "<<fname);
	    IO_.DumpFile(fname.c_str(), ss);
            ++last_checkpoint_;
        }
    }

    Error_t return_val(ErrorEvent::Success);
    if ( abs(A/A0_-1) > params_->error_factor  || abs(V/V0_-1) > params_->error_factor )
        return_val = ErrorEvent::AccuracyError;

    return return_val;
}
