template<typename VecContainer>
VesInteraction<VecContainer>::VesInteraction(int num_threads) :
    num_threads_(num_threads),
    each_thread_np_(new size_t[num_threads_]),
    each_thread_idx_(new size_t[num_threads_ + 1]),
    np_(0),
    containers_capacity_(0),
    all_pos_(NULL),
    all_den_(NULL),
    all_pot_(NULL)
{
    each_thread_idx_[0] = 0;
}

template<typename VecContainer>
VesInteraction<VecContainer>::~VesInteraction()
{
    delete[] each_thread_np_;
    delete[] each_thread_idx_;

    delete[] all_pos_;
    delete[] all_den_;
    delete[] all_pot_;
}    

template<typename VecContainer>
void VesInteraction<VecContainer>::operator()(const VecContainer &position, 
    VecContainer &density, VecContainer &potential)
{
    size_t np(position.getNumSubs() * position.getStride());
    size_t n_cpy(position.size());
    size_t idx(getCpyDestIdx(np));
    
    if(typeid(value_type) == typeid(fmm_value_type))
    {
        position.getDevice().Memcpy(all_pos_ + idx, position.begin(),
            n_cpy * sizeof(value_type), MemcpyDeviceToHost);
        
        density.getDevice().Memcpy(all_den_ + idx, density.begin(),
            n_cpy * sizeof(value_type), MemcpyDeviceToHost);
    }
    else
    {
        value_type* buffer(new value_type[n_cpy]);
        
        position.getDevice().Memcpy(buffer, position.begin(),
            n_cpy * sizeof(value_type), MemcpyDeviceToHost);
        
        for(size_t ii=0; ii<n_cpy; ++ii)
            *(all_pos_ + idx + ii) = static_cast<fmm_value_type>(buffer[ii]);
        
        position.getDevice().Memcpy(buffer, density.begin(),
            n_cpy * sizeof(value_type), MemcpyDeviceToHost);
        
        for(size_t ii=0; ii<n_cpy; ++ii)
            *(all_den_ + idx + ii) = static_cast<fmm_value_type>(buffer[ii]);

        delete[] buffer;
    }
    
    if ( omp_get_thread_num() == 0 )
        fmmInteraction();
    
    if(typeid(value_type) == typeid(fmm_value_type))
    {
        potential.getDevice().Memcpy(potential.begin(), all_pot_ + idx,
            n_cpy * sizeof(value_type), MemcpyHostToDevice);
    }
    else
    {
        value_type* buffer(new value_type[n_cpy]);
        
        for(size_t ii=0; ii<n_cpy; ++ii)
            buffer[ii] = static_cast<fmm_value_type>(*(all_pot_ + idx + ii));
   
        potential.getDevice().Memcpy(potential.begin(), buffer,
            n_cpy * sizeof(value_type), MemcpyHostToDevice);
        
        delete[] buffer;
    }
}

template<typename VecContainer>
size_t VesInteraction<VecContainer>::getCpyDestIdx(size_t this_thread_np)
{   
    int threadNum = omp_get_thread_num();
    each_thread_np_[threadNum] = this_thread_np;

#pragma omp barrier
    {
        if ( threadNum == 0 )
        {
            np_ = 0;
            for(int ii=1; ii<=num_threads_; ++ii)
            {
                np_ += each_thread_np_[ii-1];
                each_thread_idx_[ii] = each_thread_idx_[ii-1] + 
                    VecContainer::getTheDim() * each_thread_np_[ii-1]; 
            }
        }
    }

#pragma omp barrier
    {
        checkContainersSize();
    }

    return(each_thread_idx_[threadNum]);
}

template<typename VecContainer>
void VesInteraction<VecContainer>::checkContainersSize()
{
#pragma omp master
    {
        size_t new_capacity(each_thread_idx_[num_threads_]); 

        if ( containers_capacity_ < new_capacity )
        {

            delete[] all_pos_;
            all_pos_ = new fmm_value_type[new_capacity];

            delete[] all_den_;
            all_den_ = new fmm_value_type[new_capacity];
            
            delete[] all_pot_;
            all_pot_ = new fmm_value_type[new_capacity];
            
            containers_capacity_ = new_capacity;
        }
    }
#pragma omp barrier
}


template<typename VecContainer>
void VesInteraction<VecContainer>::fmmInteraction() const
{
    cout<<"FMM"<<endl;
}
//        VecContainer wrk1, wrk2, wrk3;
//         wrk1.replicate(position);
//         wrk2.replicate(position);
//         wrk3.replicate(position);
        
//         size_t np(position.getNumSubs() * position.getStride());
        
//         position.getDevice().Transpose(position.begin(), np, 3, wrk1.begin());
//         position.getDevice().Transpose( density.begin(), np, 3, wrk2.begin());
        
//         position.getDevice().DirectStokes(wrk1.begin(), wrk2.begin(), NULL, 
//             np, 1, wrk1.begin(), 0, np, wrk3.begin());  
        
//         position.getDevice().Transpose(wrk3.begin(), 3, np, potential.begin());
 
