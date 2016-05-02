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
    const typename EvolveSurface::Sca_t::device_type& device=EvolveSurface::Sca_t::getDevice();
    size_t N_ves=state->S_->getPosition().getNumSubs();

    area_new.replicate(state->S_->getPosition());
    vol_new .replicate(state->S_->getPosition());
    state->S_->area  (area_new);
    state->S_->volume( vol_new);

    if(A0_ < 0){ // Initialize area0_, vol0_
        area0_.replicate(state->S_->getPosition());
        vol0_ .replicate(state->S_->getPosition());
        device.Memcpy(area0_.begin(), area_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
        device.Memcpy( vol0_.begin(),  vol_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
        A0_=(device.MaxAbs(area0_.begin(), N_ves));
        V0_=(device.MaxAbs( vol0_.begin(), N_ves));
    }

    device.axpy(-1.0, area0_.begin(), area_new.begin(), N_ves, area_new.begin());
    device.axpy(-1.0,  vol0_.begin(),  vol_new.begin(), N_ves,  vol_new.begin());
    value_type DA(device.MaxAbs(area_new.begin(), N_ves));
    value_type DV(device.MaxAbs( vol_new.begin(), N_ves));

#pragma omp critical (monitor)
    {
        INFO(emph<<"Monitor: thread = "<<omp_get_thread_num()<<"/"<<omp_get_num_threads()
             <<", progress = "<<static_cast<int>(100*t/params_->time_horizon)<<"%"
             <<", t = "<<SCI_PRINT_FRMT<<t
             <<", dt = "<<SCI_PRINT_FRMT<<dt
             <<", area error = "<<SCI_PRINT_FRMT<<(DA/A0_)
             <<", volume error = "<<SCI_PRINT_FRMT<<(DV/V0_)<<emph);


        int checkpoint_index(checkpoint_stride_ <= 0 ? last_checkpoint_+1 : t/checkpoint_stride_);

        if ( checkpoint_flag_ && checkpoint_index > last_checkpoint_ )
        {
            ++time_idx_;
            std::string fname(params_->checkpoint_file_name);
            char suffix[6];
            sprintf(suffix, "%06d", time_idx_);
            d_["time_idx"] = std::string(suffix);
            expand_template(&fname, d_);

            //This order of packing is used when loading checkpoints
            //in ves3d_simulation file
            std::stringstream ss;
            ss<<std::scientific<<std::setprecision(16);
            params_->pack(ss, Streamable::ASCII);
            state->pack(ss, Streamable::ASCII);

            INFO("Writing data to file "<<fname);
            IO_.DumpFile(fname.c_str(), ss);
            ++last_checkpoint_;

#if HAVE_PVFMM
            if(params_->write_vtk){
                std::string vtkfbase("snapshot_");
                vtkfbase += suffix;
                INFO("Writing VTK file");

                const typename EvolveSurface::Sur_t *S_up;
                CHK( state->getSurfaceUp(S_up));
                WriteVTK(*S_up,vtkfbase.c_str());
            }
#endif // HAVE_PVFMM
        }
    }

    Error_t return_val(ErrorEvent::Success);
    if ( (DA/A0_) > params_->error_factor  || (DV/V0_) > params_->error_factor )
        return_val = ErrorEvent::AccuracyError;

    return return_val;
}
