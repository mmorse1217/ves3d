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

class DataIO
{
  public:
    explicit DataIO(string file_name = "", size_t buffer_size = 0, 
        int resize_factor = 2);    
    ~DataIO();
    
    template<typename Container>
    bool ReadData(const string &file_name, Container &data) const;
    
    template<typename Container>
    bool WriteData(const string &file_name, const Container &data, 
        ios_base::openmode mode = ios::out) const;
    
    template<typename Container>
    bool Append(const Container &data) const; 

    // Basic type IO
    template<typename T>
    bool ReadData(const string &file_name, size_t size, T* data) const;
    
    template<typename T>
    bool WriteData(const string &file_name, size_t size, const T* data, 
        ios_base::openmode mode = ios::out) const;
    
    template<typename T>
    bool Append(const T *data, size_t length) const; 

  private:
    string out_file_;
    mutable size_t out_size_;
    mutable size_t out_used_;
    mutable char* out_buffer_; 
    int resize_factor_;
   
    bool ResizeOutBuffer(size_t buffer_size) const;
    bool FlushBuffer() const;
};

#include "DataIO_templates.cc"

#endif

    

