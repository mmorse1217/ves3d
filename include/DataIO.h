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
#include <string>
#include <cassert>
#include "Logger.h"
#include "enums.h"
/**
 * A simple data I/O class tailored for this project. It reads data
 * from and writes data to file. 
 */
using namespace std;

template<typename T, typename Device>
class DataIO
{
  private:
    mutable size_t out_size_;
    mutable size_t out_used_;
    string out_file_name_;
    mutable T* out_buffer_; 
    int resize_factor_;

    const Device &device_;
  public:
    DataIO(const Device &device_in, string file_name_in = "", 
        size_t buffer_size_in = 0);
    ~DataIO();

    // Basic type IO
    bool ReadData(const char* file_name_in, size_t data_size_in, 
        T* data_array_out) const;
    bool WriteData(const char *file_name_in, size_t data_size_in, 
        const T* data_array_in, ios_base::openmode mode_in = ios::out) const;
    bool Append(const T *x_in, size_t length) const; 
    bool ResizeOutBuffer(size_t size_in) const;
    bool FlushBuffer() const;

    // General container IO
    template<typename Container>
    bool ReadData(const char* file_name_in, Container &data) const;
    
    template<typename Container>
    bool WriteData(const char *file_name_in, const Container &data, 
        ios_base::openmode mode_in = ios::out) const;
    
    template<typename Container>
    bool Append(const Container &data) const; 
};

#include "DataIO.cc"

#endif

    

