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
#include <memory>
#include <vector>

#include "Logger.h"
#include "Error.h"
#include "ves3d_common.h"

/**
 * A simple data I/O class tailored for this project. It reads data
 * from and writes data to file.
 *
 * @note: Make sure flush the buffer with the type you'd like. If this
 * is deconstructed with non-empty buffer, the buffer is dumped as
 * bindary.
 */
class DataIO
{
  public:
    enum IOFormat {BIN,ASCII};

    explicit DataIO(std::string file_name = "", IOFormat frmt = ASCII,
        size_t buffer_size = 0, int resize_factor = 2);
    ~DataIO();

    template<typename Container>
    bool ReadData(const std::string &file_name, Container &data,
        IOFormat frmt = ASCII, int offset = 0 , int length = -1) const;

    template<typename T>
    bool ReadData(const std::string &file_name, std::vector<T> &data,
        IOFormat frmt = ASCII);

    template<typename Container>
    bool WriteData(const std::string &file_name, const Container &data,
        IOFormat frmt = ASCII, std::ios_base::openmode mode = std::ios::out) const;

    template<typename Container>
    bool Append(const Container &data) const;

    bool ResizeOutBuffer(size_t buffer_size) const;

    template<typename T>
    bool FlushBuffer() const;

    static Error_t SlurpFile(const char* fname, std::ostream &content);
    static Error_t DumpFile(const char* fname, std::ostream &content);

  private:
    // Basic type IO
    // IOFormat default is differnet from public methods b/c of legacy
    template<typename T>
    bool ReadData(const std::string &file_name, size_t size,
        T* data, IOFormat frmt = BIN) const;

    template<typename T>
    bool WriteData(const std::string &file_name, size_t size, const T* data,
        IOFormat frmt = BIN,  std::ios_base::openmode mode = std::ios::out) const;

    bool FlushBufferBin() const;
    bool ReadBin(const std::string &file_name, size_t size, char* data) const;
    bool WriteBin(const std::string &file_name, size_t size, const char* data,
        std::ios_base::openmode mode = std::ios::out) const;

    std::string out_file_;
    IOFormat frmt_;
    mutable size_t out_size_;
    mutable size_t out_used_;
    mutable char* out_buffer_;
    int resize_factor_;
};

std::string FullPath(const std::string fname);
std::string FullPath(const char* fname);

#include "DataIO_templates.cc"

#endif
