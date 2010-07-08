template<>
bool DataIO::ReadData(const string &file_name, size_t size, char* data) const
{
    ifstream data_file(file_name.c_str(), ios::binary | ios::in);
    
    if(!data_file.good())
        CERR("\n Could not read the data from the file."
            <<"\n\n File name : "<<file_name<<".\n",endl, exit(1));
    
    data_file.read(data, size); 
    data_file.close();
    return(true);
}


template<typename T>
bool DataIO::ReadData(const string &file_name, size_t size, T* data) const
{
    ifstream data_file(file_name.c_str(), ios::in);
    
    if(!data_file.good())
        CERR("\n Could not read the data from the file."
            <<"\n\n File name : "<<file_name<<".\n",endl, exit(1));
    
    size_t idx=0;
    while (idx<size )
        data_file>>data[idx++];
    
    data_file.close();
    return(true);
}

template<>
bool DataIO::WriteData(const string &file_name, size_t size, const char* data, 
    ios_base::openmode mode) const
{
    ofstream data_file(file_name.c_str(), ios::binary | mode);
    if(!data_file)
        CERR(" Could not write the data to the file." 
            <<"\n\n File name : "<< file_name, endl, exit(1));
    data_file.write(data, size);    
    data_file.close();
    return(true);
}

template<typename T>
bool DataIO::WriteData(const string &file_name, size_t size, const T* data, 
    ios_base::openmode mode) const
{
    ofstream data_file(file_name.c_str(), mode);
    if(!data_file)
        CERR(" Could not write the data to the file." 
            <<"\n\n File name : "<< file_name, endl, exit(1));
    
    size_t idx=0;
    while  (idx <size )
        data_file<<data[idx++]<<endl;
    
    data_file.close();
    return(true);
}

template<typename T>
bool DataIO::Append(const T* data, size_t length) const
{
    COUTDEBUG("\n  DataIO::Append():"
        <<"\n              Size      = "<<length
        <<"\n              Total     = "<<out_size_
        <<"\n              Used      = "<<out_used_
        <<"\n              Available = "<<out_size_-out_used_<<endl);

    length *= sizeof(T);
    
    if(length > out_size_)
        ResizeOutBuffer(length);
    
    if(length > (out_size_ - out_used_))
        FlushBuffer();
    
    memcpy(out_buffer_ + out_used_, data, length);
    out_used_ +=length;
    
    return(true);
}


template<typename Container>
bool DataIO::ReadData(const string &file_name, Container &data) const
{
    size_t size = data.size() * sizeof(typename Container::value_type);

    if(Container::getDeviceType() == CPU)
        return(this->ReadData(file_name, size, (char*) data.begin()));
    else
    {
        char* buffer =new char[size];
        bool ret = ReadData(file_name, size, buffer);
        
        data.getDevice().Memcpy(data.begin(), buffer, size
            , MemcpyHostToDevice);
        
        delete[] buffer;
        return  ret;
    }
}

template<typename Container>
bool DataIO::WriteData(const string &file_name, const Container &data, 
    ios_base::openmode mode) const
{
    return(this->WriteData(file_name, data.size(), data.begin(), mode));
}
    
template<typename Container>
bool DataIO::Append(const Container &data) const
{
    return(this->Append(data.begin(), data.size()));
}
