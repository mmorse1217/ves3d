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

DataIO::~DataIO()
{
    if(out_used_ > 0){
        FlushBufferBin();
        CERR_LOC("DataIO object deconstructed with non-empty buffer"
            ", dumping the content as binary","",NULL);
    }

    free(out_buffer_);
}

bool DataIO::ReadBin(const std::string &file_name, size_t size, char* data) const
{
    COUTDEBUG("Reading binary file: "<<file_name);

    std::ifstream data_file(file_name.c_str(), std::ios::binary | std::ios::in);

    if(!data_file.good())
        CERR_LOC("Cannot open the file for reading: "<<file_name,"", exit(1));

    data_file.read(data, size);
    data_file.close();
    return(true);
}

bool DataIO::WriteBin(const std::string &file_name, size_t size, const char* data,
    std::ios_base::openmode mode) const
{
#pragma omp ordered //critical writeData
    {
        COUTDEBUG("Writing to binary file '"<<file_name<<"'"<<", size = "<<size);

        std::ofstream data_file(file_name.c_str(), std::ios::binary | mode);
        if(!data_file)
            CERR_LOC("Cannot open file for writing: "<<file_name, "", exit(1));

        data_file.write(data, size);
        data_file.close();
    }
    return(true);
}

bool DataIO::FlushBufferBin() const
{
    bool res(true);
    if(out_buffer_ !=0 && out_used_ > 0)
    {
        COUTDEBUG("Flush buffer (char) to: "<<out_file_<<", length="<<out_used_);
        res = this->WriteBin(out_file_, out_used_, out_buffer_, std::ios::app);
        out_used_ = 0;
    }

    return(res);
}

bool DataIO::ResizeOutBuffer(size_t buffer_size) const
{
    this->FlushBuffer<char>();
    free(out_buffer_);
    out_buffer_ = (char*) malloc(buffer_size);
    out_size_ = buffer_size;
    return(true);
}

Error_t DataIO::SlurpFile(const char* fname, std::ostream &content)
{
    std::ifstream fh(fname, std::ios::in);

    if(!fh)
	CERR_LOC("Cannot open file for reading: "<<fname, "", exit(1));

    content<<fh.rdbuf();
    fh.close();

    return ErrorEvent::Success;
}

Error_t DataIO::DumpFile(const char* fname, std::ostream &content)
{
    std::ofstream fh(fname, std::ios::out);

    COUT(fname);
    if(!fh)
	CERR_LOC("Cannot open file for writing: "<<fname, "", exit(1));

    fh<<content.rdbuf();
    fh.close();

    return ErrorEvent::Success;
}

std::string FullPath(const std::string fname){
    std::string base(VES3D_PATH);
    base += "/" + fname;
    return base;
}

std::string FullPath(const char* fname){
    std::string str(fname);
    return FullPath(str);
}
