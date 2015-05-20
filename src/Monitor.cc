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

        area_.replicate(state->S_->getPosition());
        vol_ .replicate(state->S_->getPosition());
        device.Memcpy(area_.begin(), area_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
        device.Memcpy( vol_.begin(),  vol_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
    }

    device.axpy(-1.0, area_new.begin(), area_.begin(), N_ves, area_.begin());
    device.axpy(-1.0,  vol_new.begin(),  vol_.begin(), N_ves,  vol_.begin());
    value_type dA(device.MaxAbs(area_.begin(), N_ves));
    value_type dV(device.MaxAbs( vol_.begin(), N_ves));

    device.axpy(-1.0, area_new.begin(), area0_.begin(), N_ves, area_.begin());
    device.axpy(-1.0,  vol_new.begin(),  vol0_.begin(), N_ves,  vol_.begin());
    value_type DA(device.MaxAbs(area_.begin(), N_ves));
    value_type DV(device.MaxAbs( vol_.begin(), N_ves));

    { // Copy area_new to area_, vol_new to vol_
      area_.replicate(state->S_->getPosition());
      vol_ .replicate(state->S_->getPosition());
      device.Memcpy(area_.begin(), area_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
      device.Memcpy( vol_.begin(),  vol_new.begin(), N_ves*sizeof(value_type), device.MemcpyDeviceToDevice);
    }
    value_type A(device.MaxAbs(area_.begin(), N_ves));
    value_type V(device.MaxAbs(vol_ .begin(), N_ves));

    if(A0_ > 0){ // Set new dt
       value_type timestep_order=1;
       value_type time_horizon=params_->time_horizon;

       value_type beta_A;
       { // compute beta_A
          //beta_A = (A/dA) * (dt/(time_horizon-t)) * (params_->error_factor - DA/A);
          beta_A = (A/dA) * (dt/time_horizon) * params_->error_factor; // This is more robust
          beta_A = std::pow<value_type>(beta_A, timestep_order);
       }

       value_type beta_V;
       { // compute beta_V
          //beta_V = (V/dV) * (dt/(time_horizon-t)) * (params_->error_factor - DV/V);
          beta_V = (V/dV) * (dt/time_horizon) * params_->error_factor; // This is more robust
          beta_V = std::pow<value_type>(beta_V, timestep_order);
       }

       value_type beta=std::min(beta_A, beta_V);
       beta=std::min(beta,1.5);
       beta=std::max(beta,0.2);

       value_type beta_scale=sqrt(0.9);
       dt=std::pow(beta_scale,timestep_order) * beta * dt;
    }else{ // Initialize A0_, V0_
        A0_=(device.MaxAbs(area0_.begin(), N_ves));
        V0_=(device.MaxAbs( vol0_.begin(), N_ves));
    }

#pragma omp critical (monitor)
    {

        INFO("Monitor: thread = "<<omp_get_thread_num()<<"/"<<omp_get_num_threads()
            <<", t = "<<SCI_PRINT_FRMT<<t
            <<", dt = "<<SCI_PRINT_FRMT<<dt
            <<", area error = "<<SCI_PRINT_FRMT<<(DA/A)
            <<", volume error = "<<SCI_PRINT_FRMT<<(DV/V));

        int checkpoint_index(checkpoint_stride_ <= 0 ? last_checkpoint_+1 : t/checkpoint_stride_);

        if ( checkpoint_flag_ && checkpoint_index > last_checkpoint_ )
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
    if ( (DA/A) > params_->error_factor  || (DV/V) > params_->error_factor )
        return_val = ErrorEvent::AccuracyError;

    return return_val;
}
