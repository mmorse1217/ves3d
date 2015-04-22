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
Array<T, DT, DEVICE>::Array(std::istream &is, Format format) :
    size_(0), capacity_(0), data_(NULL)
{
    //to count the number of calls
    PROFILESTART();
    Array<T, DT, DEVICE>::unpack(is, format);
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

template<typename T, typename DT, const DT &DEVICE>
Error_t Array<T, DT, DEVICE>::pack(std::ostream &os, Format format) const
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");

    const T *buffer(begin());

    if (DT::IsHost()){
	T* bb = new T[size()];
	buffer = bb;
	DEVICE.Memcpy(
	    bb,
	    begin(),
	    mem_size(),
	    DT::MemcpyDeviceToHost);
    }

    os<<"ARRAY\n";
    os<<"name: "<<name_<<"\n";
    os<<"size: "<<size()<<"\n";
    os<<"data: ";
    pack_array(os, format, buffer, size());
    os<<"\n/ARRAY\n";

    if (!DT::IsHost())
	delete[] buffer;

    return ErrorEvent::Success;
}

template<typename T, typename DT, const DT &DEVICE>
Error_t Array<T, DT, DEVICE>::unpack(std::istream &is, Format format)
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s, key;
    is>>s;
    ASSERT(s=="ARRAY", "Bad input string (missing header).");

    size_t sz(0);
    is>>key>>name_;
    ASSERT(key=="name:", "bad key");
    is>>key>>sz;
    ASSERT(key=="size:", "bad key");
    resize(sz);
    is>>key;
    ASSERT(key=="data:", "bad key");
    unpack_array(is, format, begin(), sz);
    is>>s;
    ASSERT(s=="/ARRAY", "Bad input string (missing footer).");

    return ErrorEvent::Success;
}

//////////////////////////////////////////////////////////////////////
template<typename T, typename DT, const DT &DEVICE>
std::ostream& operator<<(std::ostream& output,
    const Array<T, DT, DEVICE> &arr)
{
    arr.pack(output, Streamable::ASCII);
    return output;
}
