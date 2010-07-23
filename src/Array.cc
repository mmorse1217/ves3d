template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Array<T, DT, DEVICE>::Array(size_t size) :
    size_(size),
    capacity_(size_),
    data_((capacity_ > 0) ? (T*) DEVICE.Malloc(capacity_ * sizeof(T)) : NULL)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Array<T, DT, DEVICE>::~Array()
{
    DEVICE.Free(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
enum DeviceType Array<T, DT, DEVICE>::getDeviceType()
{
    return(DT);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
const Device<DT>& Array<T, DT, DEVICE>::getDevice()
{
    return(DEVICE);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Array<T, DT, DEVICE>::size() const
{
    return(size_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Array<T, DT, DEVICE>::resize(size_t new_size)
{
    ///Note that new_size = 0 frees all the memory (no recovery).
    if ( new_size > capacity_ )
    {
        T *data_new((T*) DEVICE.Malloc(new_size * sizeof(T)));
        
        if ( data_ != NULL ) 
        {
            DEVICE.Memcpy(data_new, data_, this->size() * sizeof(T), 
                MemcpyDeviceToDevice);
            DEVICE.Free(data_);
        }
        data_ = data_new;
        capacity_ = new_size;
    } 
    else if ( new_size == 0 && data_ != NULL )
    {
        DEVICE.Free(data_);
        data_ = NULL;
        capacity_ = 0;
    }

    size_ = new_size;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Array<T, DT, DEVICE>::iterator Array<T, DT, DEVICE>::begin()
{
    return(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Array<T, DT, DEVICE>::const_iterator Array<T, DT, DEVICE>::begin() const
{
    return(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Array<T, DT, DEVICE>::iterator Array<T, DT, DEVICE>::end()
{
    return(data_ + this->size());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Array<T, DT, DEVICE>::const_iterator Array<T, DT, DEVICE>::end() const
{
    return(data_ + this->size());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
bool Array<T, DT, DEVICE>::isNan() const
{
    return(DEVICE.isNan(data_, this->size()));
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Array<T, DT, DEVICE>::fillRand()
{
    return(DEVICE.fillRand(data_, this->size()));
}
