#include "DataIO.h"

DataIO::DataIO(string file_name, size_t buffer_size, 
    int resize_factor) :
    out_file_(file_name),
    out_size_(buffer_size),
    out_used_(0),
    out_buffer_((char*) malloc(out_size_)),
    resize_factor_(resize_factor)
{}
    
DataIO::~DataIO()
{
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
    COUTDEBUG("\n  DataIO::FlushBuffer() to:\n"
        <<"              "<<out_file_ <<endl);
    
    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        res = this->WriteData(out_file_, out_used_, out_buffer_, ios::app);
        out_used_ = 0;
    }
    
    return(res);
}
    

