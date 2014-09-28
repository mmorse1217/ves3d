#include "DataIO.h"

DataIO::DataIO(std::string file_name, IOFormat frmt,
    size_t buffer_size, int resize_factor) :
    out_file_(file_name),
    frmt_(frmt),
    out_size_(buffer_size),
    out_used_(0),
    out_buffer_((char*) malloc(out_size_)),
    resize_factor_(resize_factor)
{}

template<>
bool DataIO::ReadData<char>(const std::string &file_name, size_t size,
    char* data, IOFormat frmt) const
{
    COUTDEBUG("Reading binary file: "<<file_name);

    ASSERT(frmt==BIN,"char IO is only supported for BIN format");
    std::ifstream data_file(file_name.c_str(), std::ios::binary | std::ios::in);

    if(!data_file.good())
        CERR_LOC("Could not read the data from the file: "<<file_name,"", exit(1));

    data_file.read(data, size);
    data_file.close();
    return(true);
}

template<>
bool DataIO::WriteData<char>(const std::string &file_name, size_t size, const char* data,
    IOFormat frmt, std::ios_base::openmode mode) const
{
#pragma omp ordered //critical writeData
    {
        COUTDEBUG("Writing to binary file '"<<file_name<<"'"<<", size = "<<size);

        ASSERT(frmt==BIN,"char IO is only supported for BIN format");
        std::ofstream data_file(file_name.c_str(), std::ios::binary | mode);
        if(!data_file)
            CERR_LOC("Could not write the data to the file: "<<file_name, "", exit(1));

        data_file.write(data, size);
        data_file.close();
    }
    return(true);
}

template<>
bool DataIO::FlushBuffer<char>() const
{
    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        COUTDEBUG("Flush buffer (char) to: "<<out_file_<<", length="<<out_used_);
        res = this->WriteData(out_file_, out_used_, out_buffer_, BIN, std::ios::app);
        out_used_ = 0;
    }

    return(res);
}

DataIO::~DataIO()
{
    if(out_used_ > 0){
        FlushBuffer<char>();
        CERR_LOC("DataIO object deconstructed with non-empty buffer"
            ", dumping the content as binary","",NULL);
    }

    free(out_buffer_);
}

bool DataIO::ResizeOutBuffer(size_t buffer_size) const
{
    this->FlushBuffer<char>();
    free(out_buffer_);
    out_buffer_ = (char*) malloc(buffer_size);
    out_size_ = buffer_size;
    return(true);
}
