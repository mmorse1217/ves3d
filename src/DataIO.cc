#include "DataIO.h"

DataIO::DataIO(std::string file_name, size_t buffer_size,
    int resize_factor) :
    out_file_(file_name),
    out_size_(buffer_size),
    out_used_(0),
    out_buffer_((char*) malloc(out_size_)),
    resize_factor_(resize_factor)
{}

DataIO::~DataIO()
{
    if(out_used_ > 0)
        FlushBuffer();
    free(out_buffer_);
}

bool DataIO::ResizeOutBuffer(size_t buffer_size) const
{
    this->FlushBuffer();
    free(out_buffer_);
    out_buffer_ = (char*) malloc(buffer_size);
    out_size_ = buffer_size;
    return(true);
}

bool DataIO::FlushBuffer() const
{
    COUT("\n  DataIO::FlushBuffer() to:\n"
        <<"              "<<out_file_ <<std::endl);

    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        res = this->WriteData(out_file_, out_used_, out_buffer_, std::ios::app);
        out_used_ = 0;
    }

    return(res);
}

bool DataIO::ReadData(const std::string &file_name, size_t size, char* data) const
{
    std::ifstream data_file(file_name.c_str(), std::ios::binary | std::ios::in);

    if(!data_file.good())
        CERR("\n Could not read the data from the file."
            <<"\n\n File name : "<<file_name<<".\n",std::endl, exit(1));

    COUT("\n  Reading file:\n    "<<file_name<<std::endl<<std::endl);

    data_file.read(data, size);
    data_file.close();
    return(true);
}

bool DataIO::WriteData(const std::string &file_name, size_t size, const char* data,
    std::ios_base::openmode mode) const
{
    std::ofstream data_file(file_name.c_str(), std::ios::binary | mode);
    if(!data_file)
        CERR(" Could not write the data to the file."
            <<"\n\n File name : "<< file_name, std::endl, exit(1));
    data_file.write(data, size);
    data_file.close();
    return(true);
}
