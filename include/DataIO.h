/**
 * @file   DataIO.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Thu Feb  4 13:11:15 2010
 * 
 * @brief  The data I/O class header.
 */

#ifndef _DATAIO_H_
#define _DATAIO_H_

#include "Device.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>

/**
 * A simple data I/O class tailored for this project. It reads data
 * from and writes data to file. 
 */
using namespace std;

template<typename T,enum DeviceType DT>
class DataIO
{
  private:
    mutable size_t out_size_;
    mutable size_t out_used_;
    string out_file_name_;
    //Device<DT> *device_;
    mutable T* out_buffer_; 
    int resize_factor_;

  public:
    ///@todo device should be private, be OperatorMats needs it!!
    Device<DT> *device_;
  public:
    DataIO(Device<DT> *device_in, string file_name_in = "", 
        size_t buffer_size_in = 0);
    ~DataIO();
    
    bool ReadData(const char* file_name_in, size_t data_size_in, 
        T* data_array_out) const;
    
    bool WriteData(const char *file_name_in, size_t data_size_in, 
        const T* data_array_in, ios_base::openmode mode_in = ios::out) const;

    bool Append(const T *x_in, size_t length) const; 
    bool ResizeOutBuffer(size_t size_in) const;
    bool FlushBuffer() const;
};

#include "DataIO.cc"

#endif

    

