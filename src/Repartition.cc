template<typename T>
Repartition<T>::Repartition(GlobalRepart_t fun_ptr,
    Dealloc_t clear_context,
    int num_threads) :
    g_repart_handle_(fun_ptr),
    clear_context_(clear_context),
    num_threads_(num_threads),
    each_thread_nv_(new size_t[num_threads_]),
    each_thread_idx_(new size_t[num_threads_ + 1]),
    nv_(0),
    capacity_(0),
    all_pos_(NULL),
    all_tension_(NULL),
    posr_(NULL),
    tensionr_(NULL),
    context_(NULL),
    nvr_(0)
{
    COUTDEBUG("creating a repartion object");
    for(int ii=0; ii<num_threads_; ++ii)
        each_thread_idx_[ii] = each_thread_nv_[ii] = 0;

    if (this->g_repart_handle_ && !this->clear_context_)
      WARN("No deallocator is defined for the repartition_context."
	" This may cause memory leak.");
}

template<typename T>
Repartition<T>::~Repartition()
{
    COUTDEBUG("destroying the repartion object");
    delete[] each_thread_nv_;
    delete[] each_thread_idx_;

    delete[] all_pos_;
    delete[] all_tension_;

    if (this->context_){
        COUTDEBUG("deleting the repartion context");
        this->clear_context_(&(this->context_));
    }
}


template<typename T>
template<typename VecContainer, typename ScaContainer>
Error_t Repartition<T>::operator()(VecContainer &coord,
    ScaContainer &tension) const
{
    assert( typeid(T) == typeid(typename VecContainer::value_type) );
    assert( typeid(T) == typeid(typename ScaContainer::value_type) );
    assert( coord.getNumSubs() == tension.getNumSubs() );

    if ( g_repart_handle_ == NULL )
    {
        COUTDEBUG("No repartition handle, returning");
        return(ErrorEvent::ReferenceError);
    }

    //Getting the sizes
    size_t nv(tension.getNumSubs());
    size_t stride(tension.getStride());
    size_t idx(this->getCpyIdx(nv, stride));

#pragma omp barrier
    {
        checkContainersSize(stride);
    }

    //Copying to the host
    VecContainer::getDevice().Memcpy(all_pos_ + VecContainer::getTheDim() * idx,
        coord.begin(), coord.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    ScaContainer::getDevice().Memcpy(all_tension_ + idx, tension.begin(),
        tension.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    // call user interaction routine
#pragma omp barrier

#pragma omp master
    COUTDEBUG("repartitioning vesicle distribution with "<<nv_<<" vesicles");
    g_repart_handle_(nv_, stride, all_pos_, all_tension_, &nvr_,
        &posr_, &tensionr_, &(this->context_));

#pragma omp barrier

    int oldnv(nv);
    nv = getNvShare();
    idx = this->getCpyIdx(nv, stride);
    coord.resize(nv);
    tension.resize(nv);

    //Copying back the new values to the device(s)
    VecContainer::getDevice().Memcpy(coord.begin(), posr_ + VecContainer::getTheDim() * idx,
        coord.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    ScaContainer::getDevice().Memcpy(tension.begin(), tensionr_ + idx,
        tension.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

#pragma omp master
    {
        delete[] posr_;
        delete[] tensionr_;
    }


#pragma omp critical (reparamPrint)
    {
        COUTDEBUG("Repartitioning thread = "
            <<omp_get_thread_num()<<"/"<<omp_get_num_threads()
            <<", initial surfaces = "<<oldnv
            <<", new surfaces = "<<nv
                  );
    }

    return(ErrorEvent::Success);
}

template<typename T>
template<typename VecContainer, typename ScaContainer, typename ArrContainer>
Error_t Repartition<T>::operator()(VecContainer &coord, ScaContainer &tension, ArrContainer &contrast, 
        ArrContainer &excess_density, ArrContainer &bending_modulus, VecContainer &vel, VecContainer &fc) const
{
    assert( typeid(T) == typeid(typename VecContainer::value_type) );
    assert( typeid(T) == typeid(typename ScaContainer::value_type) );
    assert( typeid(T) == typeid(typename ArrContainer::value_type) );
    assert( coord.getNumSubs() == tension.getNumSubs() );
    assert( coord.getNumSubs() == vel.getNumSubs() );
    assert( coord.getNumSubs() == fc.getNumSubs() );

    if ( g_repart_handle_ == NULL )
    {
        COUTDEBUG("No repartition handle, returning");
        return(ErrorEvent::ReferenceError);
    }

    //Getting the sizes
    size_t nv(tension.getNumSubs());
    size_t stride(tension.getStride());
    size_t idx(this->getCpyIdx(nv, stride));

#pragma omp barrier
    {
        checkContainersSize(stride);
    }

    //Copying to the host
    VecContainer::getDevice().Memcpy(all_pos_ + VecContainer::getTheDim() * idx,
        coord.begin(), coord.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    ScaContainer::getDevice().Memcpy(all_tension_ + idx, tension.begin(),
        tension.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    ArrContainer::getDevice().Memcpy(all_contrast_ + idx/stride,
        contrast.begin(), contrast.size() * sizeof(T),
        ArrContainer::getDevice().MemcpyDeviceToHost);

    ArrContainer::getDevice().Memcpy(all_excess_density_ + idx/stride,
        excess_density.begin(), excess_density.size() * sizeof(T),
        ArrContainer::getDevice().MemcpyDeviceToHost);

    ArrContainer::getDevice().Memcpy(all_bending_modulus_ + idx/stride,
        bending_modulus.begin(), bending_modulus.size() * sizeof(T),
        ArrContainer::getDevice().MemcpyDeviceToHost);

    VecContainer::getDevice().Memcpy(all_vel_ + VecContainer::getTheDim() * idx,
        vel.begin(), vel.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    VecContainer::getDevice().Memcpy(all_fc_ + VecContainer::getTheDim() * idx,
        fc.begin(), fc.size() * sizeof(T),
        VecContainer::getDevice().MemcpyDeviceToHost);

    // call user interaction routine
#pragma omp barrier

#pragma omp master
    {
        COUTDEBUG("repartitioning vesicle distribution with "<<nv_<<" vesicles");
        g_repart_handle_(nv_, stride, 
                all_pos_, all_tension_, all_contrast_, all_excess_density_, all_bending_modulus_, all_vel_, all_fc_,
                &nvr_,
                &posr_, &tensionr_, contrastr_, excess_densityr_, bending_modulusr_, velr_, fcr_, 
                &(this->context_));
    }

#pragma omp barrier

    int oldnv(nv);
    nv = getNvShare();
    idx = this->getCpyIdx(nv, stride);
    coord.resize(nv);
    tension.resize(nv);
    contrast.resize(nv);
    excess_density.resize(nv);
    bending_modulus.resize(nv);
    vel.resize(nv);
    fc.resize(nv);

    //Copying back the new values to the device(s)
    VecContainer::getDevice().Memcpy(coord.begin(), posr_ + VecContainer::getTheDim() * idx,
        coord.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    ScaContainer::getDevice().Memcpy(tension.begin(), tensionr_ + idx,
        tension.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    ArrContainer::getDevice().Memcpy(contrast.begin(), contrastr_ + idx/stride,
        contrast.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    ArrContainer::getDevice().Memcpy(excess_density.begin(), excess_densityr_ + idx/stride,
        excess_density.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    ArrContainer::getDevice().Memcpy(bending_modulus.begin(), bending_modulusr_ + idx/stride,
        bending_modulus.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    VecContainer::getDevice().Memcpy(vel.begin(), velr_ + VecContainer::getTheDim() * idx,
        vel.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

    VecContainer::getDevice().Memcpy(fc.begin(), fcr_ + VecContainer::getTheDim() * idx,
        fc.size() * sizeof(T),
        VecContainer::getDevice().MemcpyHostToDevice);

#pragma omp master
    {
        delete[] posr_;
        delete[] tensionr_;
        delete[] contrastr_;
        delete[] excess_densityr_;
        delete[] bending_modulusr_;
        delete[] velr_;
        delete[] fcr_;
    }


#pragma omp critical (reparamPrint)
    {
        COUTDEBUG("Repartitioning thread = "
            <<omp_get_thread_num()<<"/"<<omp_get_num_threads()
            <<", initial surfaces = "<<oldnv
            <<", new surfaces = "<<nv
                  );
    }

    return(ErrorEvent::Success);
}

template<typename T>
size_t Repartition<T>::getCpyIdx(size_t this_thread_nv,
    size_t stride) const
{
    int threadNum = omp_get_thread_num();
    each_thread_nv_[threadNum] = this_thread_nv;
    int runtimeNumThreads = omp_get_num_threads();
#pragma omp barrier
    {
        if ( threadNum == 0 )
        {
            nv_ = 0;

            for(int ii=1; ii<=runtimeNumThreads; ++ii)
            {
                nv_ += each_thread_nv_[ii-1];
                each_thread_idx_[ii] = each_thread_idx_[ii-1] +
                    stride * each_thread_nv_[ii-1];
            }
        }
    }

#pragma omp barrier

    return(each_thread_idx_[threadNum]);
}

template<typename T>
size_t Repartition<T>::getNvShare() const
{
    int runtimeNumThreads = omp_get_num_threads();
    size_t nv(nvr_/runtimeNumThreads);

    if ( omp_get_thread_num() == 0 )
        nv = nvr_ - (runtimeNumThreads - 1) * nv;

    return nv;
}

template<typename T>
void Repartition<T>::checkContainersSize(size_t stride) const
{
#pragma omp master
    {
        if ( capacity_ < nv_ )
         {

            delete[] all_pos_;
            all_pos_ = new T[nv_ * VES3D_DIM * stride];

            delete[] all_tension_;
            all_tension_ = new T[nv_ * stride];

            delete[] all_contrast_;
            all_contrast_ = new T[nv_];

            delete[] all_excess_density_;
            all_excess_density_ = new T[nv_];

            delete[] all_bending_modulus_;
            all_bending_modulus_ = new T[nv_];

            delete[] all_vel_;
            all_vel_ = new T[nv_ * VES3D_DIM * stride];

            delete[] all_fc_;
            all_fc_ = new T[nv_ * VES3D_DIM * stride];

            capacity_ = nv_;
        }
    }
#pragma omp barrier
}
