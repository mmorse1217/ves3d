template<typename T>
bool DataIO::ReadData(const std::string &file_name, size_t size,
    T* data, IOFormat frmt) const
{
    if (frmt == BIN)
        return(ReadBin(file_name, size * sizeof(T),
                reinterpret_cast<char*>(data)));

    // ASCII input
    COUTDEBUG("Reading ASCII file: "<<file_name);
    std::ifstream data_file(file_name.c_str(), std::ios::in);

    if(!data_file.good())
        CERR_LOC("Could not read the data from the file: "<<file_name,"", exit(1));

    size_t idx=0;
    while (idx<size )
        data_file>>data[idx++];

    data_file.close();
    return(true);
}

template<typename T>
bool DataIO::WriteData(const std::string &file_name, size_t size, const T* data,
    IOFormat frmt, std::ios_base::openmode mode) const
{
    if (frmt == BIN)
        return(WriteBin(file_name, size * sizeof(T),
                reinterpret_cast<const char*>(data), mode));

#pragma omp ordered //critical writeData
    {
        COUTDEBUG("Writing to ASCII file '"<<file_name<<"'"<<", size = "<<size);
        std::ofstream data_file(file_name.c_str(), mode);
        if(!data_file)
            CERR_LOC("Could not write the data to the file: "<< file_name,"", exit(1));

        size_t idx=0;
        while  (idx <size )
            data_file<<data[idx++]<<std::endl;
        data_file.close();

        COUTDEBUG("Written "<<idx<<" items");
    }
    return(true);
}

template<typename Container>
bool DataIO::ReadData(const std::string &file_name, Container &data,
    IOFormat frmt, int offset, int length) const
{
    size_t size((length == -1)  ? data.size()-offset : length);

    if(Container::getDevice().IsHost())
        return(this->ReadData(file_name, size, data.begin()+offset, frmt));
    else
    {
        typename Container::value_type*
            buffer(new typename Container::value_type[size]);
        bool ret = ReadData(file_name, size, buffer, frmt);

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
    IOFormat frmt, std::ios_base::openmode mode) const
{
    return(this->WriteData(file_name, data.size(), data.begin(), frmt, mode));
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

        COUTDEBUG("Appending: size="<<length<<" (buffer info: total = "
            <<out_size_<<", used = "<<out_used_<<", available = "
            <<out_size_-out_used_<<")");

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
    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        COUTDEBUG("Flush buffer (typed) to: "<<out_file_<<", length="<<out_used_);
        res = this->WriteData(out_file_, out_used_/sizeof(T),
            reinterpret_cast<T*>(out_buffer_), frmt_, std::ios::app);
        out_used_ = 0;
    }

    return(res);
}
