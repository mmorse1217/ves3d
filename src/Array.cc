template<typename T, typename DT, const DT &DEVICE>
Array<T, DT, DEVICE>::Array(size_t size) :
    size_(size),
    capacity_(size_),
    data_((capacity_ > 0) ? (T*) DEVICE.Malloc(capacity_ * sizeof(T)) : NULL)
{
    //to count the number of calls
    PROFILESTART();
    PROFILEEND("",0);
}

template<typename T, typename DT, const DT &DEVICE>
Array<T, DT, DEVICE>::~Array()
{
    //mostly for counting the # of calls
    PROFILESTART();
    DEVICE.Free(this->data_);
    PROFILEEND("",0);
}

template<typename T, typename DT, const DT &DEVICE>
const DT& Array<T, DT, DEVICE>::getDevice()
{
    return(DEVICE);
}

template<typename T, typename DT, const DT &DEVICE>
size_t Array<T, DT, DEVICE>::size() const
{
    return(this->size_);
}

template<typename T, typename DT, const DT &DEVICE>
size_t Array<T, DT, DEVICE>::mem_size() const
{
    return(this->size_ * sizeof(T));
}

template<typename T, typename DT, const DT &DEVICE>
void Array<T, DT, DEVICE>::resize(size_t new_size)
{
    PROFILESTART();
    ///Note that new_size = 0 frees all the memory (no recovery).
    if ( new_size > capacity_ )
    {
        T *data_new((T*) DEVICE.Malloc(new_size * sizeof(T)));

        if ( data_ != NULL )
        {
            DEVICE.Memcpy(
                data_new,
                data_,
                this->mem_size(),
                DT::MemcpyDeviceToDevice);

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
    PROFILEEND("",0);
}

template<typename T, typename DT, const DT &DEVICE>
typename Array<T, DT, DEVICE>::iterator Array<T, DT, DEVICE>::begin()
{
    return(this->data_);
}

template<typename T, typename DT, const DT &DEVICE>
typename Array<T, DT, DEVICE>::const_iterator Array<T, DT, DEVICE>::begin() const
{
    return(this->data_);
}

template<typename T, typename DT, const DT &DEVICE>
typename Array<T, DT, DEVICE>::iterator Array<T, DT, DEVICE>::end()
{
    return(this->data_ + this->size());
}

template<typename T, typename DT, const DT &DEVICE>
typename Array<T, DT, DEVICE>::const_iterator Array<T, DT, DEVICE>::end() const
{
    return(this->data_ + this->size());
}

//////////////////////////////////////////////////////////////////////
template<typename T, typename DT, const DT &DEVICE>
std::ostream& operator<<(std::ostream& output,
    const Array<T, DT, DEVICE> &arr)
{
    typedef typename Array<T, DT, DEVICE>::value_type E;
    E *buffer(new E[arr.size()]);
    DEVICE.Memcpy(
        buffer,
        arr.begin(),
        arr.mem_size(),
        DT::MemcpyDeviceToHost);

    output<<"[";
    for(size_t ii=0; ii<arr.size(); ++ii)
        output<<buffer[ii]<<" ";
    output<<"]";

    delete[] buffer;
    return(output);
}
