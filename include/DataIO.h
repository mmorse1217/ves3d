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
  private:
    T* buffer_;
    int buffer_size_;
    int used_buffer_;
    string file_name_;
    Device<T> &device_;

  public:
    DataIO(Device<T> &device_in, string file_name_in, int buffer_size_in);
    
    ~DataIO();
    
    bool ReadData(const char* file_name_in, int data_size_in, 
        T* data_array_out);
    
    bool WriteData(const char* file_name_in, int data_size_in, 
        const T* data_array_in, ios_base::openmode mode_in = ios::out);

    bool Append(T *x_in, int length);
    
    bool ResizeBuffer(int size_in);

    bool FlushBuffer();
};


template<typename T>
bool DataIO<T>::ReadData(const char* file_name_in, 
    int data_size_in, T* data_array_out)
{
    ifstream data_file(file_name_in, ios::in);
        
    if(!data_file.good())
    {
        cerr<<"\n Could not read the data from the file." <<endl
            <<" File name : "<<file_name_in<<"\n"<<endl;
        abort();
    }

    int idx=0;
    while (idx<data_size_in)
        data_file>>data_array_out[idx++];
    
    data_file.close();
    return(true);
}

template<typename T>
bool DataIO<T>::WriteData(const char* file_name_in, 
    int data_size_in, const T* data_array_in, ios_base::openmode mode_in)
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

    int idx=0;
    while (idx<data_size_in)
        data_file<<data_array_in[idx++]<<endl;
    
    data_file.close();
    return(true);
}

template<typename T>
DataIO<T>::DataIO(Device<T> &device_in, string file_name_in, int buffer_size_in) :
    buffer_(0),
    buffer_size_(buffer_size_in),
    used_buffer_(0),
    file_name_(file_name_in),
    device_(device_in)
{
    buffer_ = (T*) malloc(buffer_size_ * sizeof(T));
}

template<typename T>
DataIO<T>::~DataIO()
{
    FlushBuffer();
    free(buffer_);
}    

template<typename T>
bool DataIO<T>::Append(T* x_in, int length)
{
#ifndef NDEBUG
    cout<<"DataIO::Append()"<<endl;
#endif

    if(length > (buffer_size_ - used_buffer_))
    {
        FlushBuffer();
    }
    
    if(length > buffer_size_)
        ResizeBuffer(length);

    device_.Memcpy(buffer_ + used_buffer_, x_in, length, MemcpyDeviceToHost);
    used_buffer_ +=length;

    return(true);
}

template<typename T>
bool DataIO<T>::ResizeBuffer(int buffer_size_in)
{
    FlushBuffer();
    free(buffer_);
    buffer_ = (T*) malloc(buffer_size_in * sizeof(T));
    buffer_size_ = buffer_size_in;
    return(true);
}

template<typename T>
bool DataIO<T>::FlushBuffer()
{
#ifndef NDEBUG
    cout<<"DataIO::FlushBuffer()"<<endl;
#endif

    if(buffer_ !=0 && used_buffer_ > 0)
    {
        bool res = WriteData(file_name_.c_str(), used_buffer_, buffer_, ios::app);
        used_buffer_ = 0;
        return(res);
    }
    else
    {
        return(true);
    }
}

#endif

    

