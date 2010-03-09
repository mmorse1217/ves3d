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

/**
 * A simple data I/O class tailored for this project. It reads data
 * from and writes data to file. 
 */
template<typename T>
class DataIO
{
  public:
    bool ReadData(const char* file_name_in, int data_size_in, 
        T* data_array_out);
    
    bool WriteData(const char* file_name_in, int data_size_in, 
        const T* data_array_in);
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
    int data_size_in, const T* data_array_in)
{
    ofstream data_file(file_name_in, ios::out);
    assert(data_file.good());

    if(!data_file)
    {
        cerr<<" Could not write the data from the file." <<endl
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

#endif

    

