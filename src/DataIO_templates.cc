template<typename T>
bool DataIO::ReadData(const string &file_name, size_t size, T* data) const
{
    ifstream data_file(file_name.c_str(), ios::in);
    
    if(!data_file.good())
        CERR("\n  Could not read the data from the file."
            <<"\n\n    File name : "<<file_name<<".\n",endl, exit(1));

    COUT("\n  Reading file:\n    "<<file_name<<endl<<endl);
    
    size_t idx=0;
    while (idx<size )
        data_file>>data[idx++];

    data_file.close();
    return(true);
}

template<typename T>
bool DataIO::WriteData(const string &file_name, size_t size, const T* data, 
    ios_base::openmode mode) const
{
#pragma omp ordered //critical writeData
    {    
        ofstream data_file(file_name.c_str(), mode);
        if(!data_file)
            CERR("  Could not write the data to the file." 
                <<"\n\n    File name : "<< file_name, endl, exit(1));
        
        size_t idx=0;
        while  (idx <size )
            data_file<<data[idx++]<<endl;
        
        data_file.close();
    }
    return(true);
}

template<typename Container>
bool DataIO::ReadData(const string &file_name, Container &data,
    int offset, int length) const
{
    size_t size((length == -1)  ? data.size()-offset : length);
    
    if(Container::getDeviceType() == CPU)
        return(this->ReadData(file_name, size, data.begin()+offset));
    else
    {
        typename Container::value_type* 
            buffer(new typename Container::value_type[size]);
        bool ret = ReadData(file_name, size, buffer);
        
        Container::getDevice().Memcpy(data.begin() + offset, buffer, 
            size * sizeof(typename Container::value_type), MemcpyHostToDevice);
        
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
    size_t length(data.size());
    length *= sizeof(typename Container::value_type);

#pragma omp critical (IOAppend)
    {
        COUTDEBUG("\n  DataIO::Append():"
            <<"\n              Size      = "<<length
            <<"\n              Total     = "<<out_size_
            <<"\n              Used      = "<<out_used_
            <<"\n              Available = "<<out_size_-out_used_<<endl);
        
        if(length > out_size_)
            ResizeOutBuffer(length);
        
        if(length > (out_size_ - out_used_))
            FlushBuffer();
        
        Container::getDevice().Memcpy(out_buffer_ + out_used_, data.begin(), 
            length, MemcpyDeviceToHost);
        
        out_used_ +=length;
    }
    return(true);
}


