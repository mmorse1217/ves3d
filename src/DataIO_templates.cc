template<typename T>
bool DataIO::ReadData(const std::string &file_name, size_t size, T* data) const
{
    std::ifstream data_file(file_name.c_str(), std::ios::in);

    if(!data_file.good())
        CERR("\n  Could not read the data from the file."
            <<"\n\n    File name : "<<file_name<<".\n",std::endl, exit(1));

    COUTDEBUG("\n  Reading file:\n    "<<file_name<<std::endl<<std::endl);

    size_t idx=0;
    while (idx<size )
        data_file>>data[idx++];

    data_file.close();
    return(true);
}

template<typename T>
bool DataIO::WriteData(const std::string &file_name, size_t size, const T* data,
    std::ios_base::openmode mode) const
{
    COUT("\n  DataIO::"<<__FUNCTION__<<":\n"
        <<"              size : "<<size<<std::endl);

#pragma omp ordered //critical writeData
    {
        std::ofstream data_file(file_name.c_str(), mode);
        if(!data_file)
            CERR("  Could not write the data to the file."
                <<"\n\n    File name : "<< file_name, std::endl, exit(1));

        size_t idx=0;
        while  (idx <size )
            data_file<<data[idx++]<<std::endl;
        data_file.close();

        COUT("           written : "<<idx<<std::endl);

    }
    return(true);
}

template<typename Container>
bool DataIO::ReadData(const std::string &file_name, Container &data,
    int offset, int length) const
{
    size_t size((length == -1)  ? data.size()-offset : length);

    if(Container::getDevice().IsHost())
        return(this->ReadData(file_name, size, data.begin()+offset));
    else
    {
        typename Container::value_type*
            buffer(new typename Container::value_type[size]);
        bool ret = ReadData(file_name, size, buffer);

        Container::getDevice().Memcpy(
            data.begin() + offset,
            buffer,
            size * sizeof(typename Container::value_type),
            Container::device_type::MemcpyHostToDevice);

        delete[] buffer;
        return  ret;
    }
}

template<typename Container>
bool DataIO::WriteData(const std::string &file_name, const Container &data,
    std::ios_base::openmode mode) const
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
        if(length > out_size_)
            ResizeOutBuffer(length);

        if(length > (out_size_ - out_used_))
            this->FlushBuffer<typename Container::value_type>();


        COUT("\n  DataIO::Append():"
            <<"\n              Size      = "<<length
            <<"\n              Total     = "<<out_size_
            <<"\n              Used      = "<<out_used_
            <<"\n              Available = "<<out_size_-out_used_
            <<std::endl);

        Container::getDevice().Memcpy(out_buffer_ + out_used_,
            data.begin(),
            length,
            Container::getDevice().MemcpyDeviceToHost);

        out_used_ +=length;
    }
    return(true);
}

template<typename T>
bool DataIO::FlushBuffer() const
{
    COUT("\n  DataIO::FlushBuffer() [typed] to:\n"
        <<"              "<<out_file_ <<std::endl);

    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        res = this->WriteData(out_file_, out_used_ /sizeof(T),
            reinterpret_cast<T*>(out_buffer_), std::ios::app);
        out_used_ = 0;
    }

    return(res);
}
