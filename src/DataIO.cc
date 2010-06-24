template<typename T, enum DeviceType DT>
DataIO<T,DT>::DataIO(const Device<DT> &device_in, 
    string file_name_in, size_t buffer_size_in) :
    out_size_(buffer_size_in),
    out_used_(0),
    out_file_name_(file_name_in),
    device_(device_in),
    out_buffer_((T*) malloc(out_size_ * sizeof(T)))
{}

template<typename T, enum DeviceType DT>
DataIO<T,DT>::~DataIO()
{
    FlushBuffer();
    free(out_buffer_);
}    

template<typename T, enum DeviceType DT>
bool DataIO<T,DT>::ReadData(const char* file_name_in, 
    size_t data_size_in, T* data_array_out) const
{
    ifstream data_file(file_name_in, ios::in);
    
    T* in_buffer = (T*) malloc(data_size_in * sizeof(T));
    if(!data_file.good())
    {
        cerr<<"\n Could not read the data from the file." <<endl
            <<" File name : "<<file_name_in<<".\n"<<endl;
        exit(1);
    }

    size_t idx=0;
    while (idx<data_size_in)
        data_file>>in_buffer[idx++];
    
    data_file.close();
    device_.Memcpy(data_array_out, in_buffer, 
        data_size_in * sizeof(T), MemcpyHostToDevice);

    free(in_buffer);
    return(true);
}

template<typename T, enum DeviceType DT>
bool DataIO<T,DT>::WriteData(const char* file_name_in, 
    size_t data_size_in, const T* data_array_in, 
    ios_base::openmode mode_in) const
{
    ofstream data_file(file_name_in, mode_in);
    assert(data_file.good());

    if(!data_file)
    {
        cerr<<" Could not write the data to the file." <<endl
            <<" File name : "<< file_name_in <<endl;
        data_file.close();
        exit(1);
    }

    size_t idx=0;
    while (idx<data_size_in)
        data_file<<data_array_in[idx++]<<endl;
    
    data_file.close();
    return(true);
}


template<typename T, enum DeviceType DT>
bool DataIO<T,DT>::Append(const T* x_in, size_t length) const
{
#ifndef NDEBUG
    COUT("DataIO::Append(): Size = "<<length<<", Available = "<<out_size_);
#endif

    if(length > out_size_)
        ResizeOutBuffer(length);

    if(length > (out_size_ - out_used_))
        FlushBuffer();
    
    device_.Memcpy(out_buffer_ + out_used_, x_in, 
        length * sizeof(T), MemcpyDeviceToHost);
    out_used_ +=length;
    
    return(true);
}

template<typename T, enum DeviceType DT>
bool DataIO<T,DT>::ResizeOutBuffer(size_t buffer_size_in) const
{
    FlushBuffer();
    free(out_buffer_);
    out_buffer_ = (T*) malloc(buffer_size_in * sizeof(T));
    out_size_ = buffer_size_in;
    return(true);
}

template<typename T, enum DeviceType DT>
bool DataIO<T,DT>::FlushBuffer() const
{
#ifndef NDEBUG
    COUT("DataIO::FlushBuffer()");
#endif

    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        res = WriteData(out_file_name_.c_str(), 
            out_used_, out_buffer_, ios::app);
        out_used_ = 0;
    }

    return(res);
}
