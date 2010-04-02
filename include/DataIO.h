/**
 * @file   DataIO.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Thu Feb  4 13:11:15 2010
 * 
 * @brief  The data I/O class header.
 */

#ifndef _DATAIO_H_
#define _DATAIO_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "Device.h"
#include <string>
#include <cassert>
/**
 * A simple data I/O class tailored for this project. It reads data
 * from and writes data to file. 
 */

using namespace std;

template<typename T>
class DataIO
{
  public:
    T* out_buffer_; 
    size_t out_size_;
    size_t out_used_;
    string out_file_name_;
    Device<T> &device_;

    bool WriteData(const char* file_name_in, size_t data_size_in, 
        const T* data_array_in, ios_base::openmode mode_in = ios::out);

  public:
    DataIO(Device<T> &device_in, string file_name_in, size_t buffer_size_in);
    
    ~DataIO();
    
    bool ReadData(const char* file_name_in, size_t data_size_in, 
        T* data_array_out);
    
    bool Append(T *x_in, size_t length);
    
    bool ResizeOutBuffer(size_t size_in);

    bool FlushBuffer();
};


template<typename T>
bool DataIO<T>::ReadData(const char* file_name_in, 
    size_t data_size_in, T* data_array_out)
{
    ifstream data_file(file_name_in, ios::in);
    
    T* in_buffer = (T*) malloc(data_size_in * sizeof(T));
    if(!data_file.good())
    {
        cerr<<"\n Could not read the data from the file." <<endl
            <<" File name : "<<file_name_in<<".\n"<<endl;
        abort();
    }

    size_t idx=0;
    while (idx<data_size_in)
        data_file>>in_buffer[idx++];
    
    data_file.close();
    device_.Memcpy(data_array_out, in_buffer, data_size_in, MemcpyHostToDevice);

    free(in_buffer);
    return(true);
}

template<typename T>
bool DataIO<T>::WriteData(const char* file_name_in, 
    size_t data_size_in, const T* data_array_in, ios_base::openmode mode_in)
{
    ofstream data_file(file_name_in, mode_in);
    assert(data_file.good());

    if(!data_file)
    {
        cerr<<" Could not write the data to the file." <<endl
            <<" File name : "<< file_name_in <<endl;
        data_file.close();
        abort();
    }

    size_t idx=0;
    while (idx<data_size_in)
        data_file<<data_array_in[idx++]<<endl;
    
    data_file.close();
    return(true);
}

template<typename T>
DataIO<T>::DataIO(Device<T> &device_in, string file_name_in, size_t buffer_size_in) :
    out_buffer_(0),
    out_size_(buffer_size_in),
    out_used_(0),
    out_file_name_(file_name_in),
    device_(device_in)
{
    out_buffer_ = (T*) malloc(out_size_ * sizeof(T));
}

template<typename T>
DataIO<T>::~DataIO()
{
    FlushBuffer();
    free(out_buffer_);
}    

template<typename T>
bool DataIO<T>::Append(T* x_in, size_t length)
{
#ifndef NDEBUG
    cout<<"DataIO::Append()"<<endl;
#endif

    if(length > out_size_)
        ResizeOutBuffer(length);

    if(length > (out_size_ - out_used_))
        FlushBuffer();
    
    device_.Memcpy(out_buffer_ + out_used_, x_in, length, MemcpyDeviceToHost);
    out_used_ +=length;
    
    return(true);
}

template<typename T>
bool DataIO<T>::ResizeOutBuffer(size_t buffer_size_in)
{
    FlushBuffer();
    free(out_buffer_);
    out_buffer_ = (T*) malloc(buffer_size_in * sizeof(T));
    out_size_ = buffer_size_in;
    return(true);
}

template<typename T>
bool DataIO<T>::FlushBuffer()
{
#ifndef NDEBUG
    cout<<"DataIO::FlushBuffer()"<<endl;
#endif

    if(out_buffer_ !=0 && out_used_ > 0)
    {
        bool res = WriteData(out_file_name_.c_str(), out_used_, out_buffer_, ios::app);
        out_used_ = 0;
        return(res);
    }
    else
    {
        return(true);
    }
}

#endif

    

