/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief unit test
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#include <typeinfo> //for typeid
#include "Logger.h"
#include "DataIO.h"
#include "Array.h"
#include "Device.h"

template<typename C>
class DataIOTest{
public:
  typedef typename C::value_type T;
  typedef typename C::device_type DT;

  bool PerformAll(){
    TestWriteReadData_ASCII();
    TestWriteReadData_BIN();
    TestAppend_ASCII();
    TestAppend_BIN();

    COUT(emph<<" *** DataIO class with "<<typeid(C).name()
	 <<" container type passed ***"<<emph);

    return true;
  }

  C* get_filled_container(size_t sz=10){

    C* X(new C(sz));
    T* D = new T [sz];
    for(size_t iD(0);iD<sz;++iD)
      D[iD] = iD;

    X->getDevice().Memcpy(X->begin(), D, sz * sizeof(T), DT::MemcpyHostToDevice);
    delete[] D;
    return X;
  }

  bool TestWriteReadData_ASCII(){
    std::string fname("DataIOTest.out");
    C* X(get_filled_container());

    DataIO io;
    // write data
    ASSERT(io.WriteData(fname,*X,DataIO::ASCII), "Write to ascii file");

    // read it back
    C Y;
    Y.resize(X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::ASCII),"Read from ascii file");

    T* B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB,"Expected content");
    delete[] B;

    // appending
    ASSERT(io.WriteData(fname,*X,DataIO::ASCII,std::ios::app), "Write to ascii file");

    // reading
    Y.resize(2*X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::ASCII),"Read from ascii file");

    B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB%10,"Expected content");


    delete[] B;

    return true;
  }


  bool TestWriteReadData_BIN(){
    std::string fname("DataIOTest.out");
    C* X(get_filled_container());

    DataIO io;
    // write data
    ASSERT(io.WriteData(fname,*X,DataIO::BIN), "Write to bin file");

    // read it back
    C Y;
    Y.resize(X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::BIN),"Read from bin file");

    T* B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB,"Expected content");
    delete[] B;

    // appending
    ASSERT(io.WriteData(fname,*X,DataIO::BIN,std::ios::app), "Write to bin file");

    // reading
    Y.resize(2*X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::BIN),"Read from bin file");

    B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB%10,"Expected content");


    delete[] B;

    return true;
  }

  bool TestAppend_ASCII(){

    std::string fname("DataIOTest.out");
    remove(fname.c_str());
    C* X(get_filled_container());

    DataIO io(fname);
    io.Append(*X);
    io.FlushBuffer<T>();

    // read it back
    C Y;
    Y.resize(X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::ASCII),"Read from ascii file");

    T* B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB,"Expected content");
    delete[] B;

    // appending
    io.Append(*X);
    io.FlushBuffer<T>();

    // reading
    Y.resize(2*X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::ASCII),"Read from ascii file");

    B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB%10,"Expected content");

    delete[] B;

    return true;
  }

  bool TestAppend_BIN(){

    std::string fname("DataIOTest.out");
    remove(fname.c_str());
    C* X(get_filled_container());

    DataIO io(fname,DataIO::BIN);
    io.Append(*X);
    io.FlushBuffer<T>();

    // read it back
    C Y;
    Y.resize(X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::BIN),"Read from ascii file");

    T* B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB,"Expected content");
    delete[] B;

    // appending
    io.Append(*X);
    io.FlushBuffer<T>();

    // reading
    Y.resize(2*X->size());
    ASSERT(io.ReadData(fname,Y,DataIO::BIN),"Read from ascii file");

    B = new T [Y.size()];
    Y.getDevice().Memcpy(B, Y.begin(), Y.size() * sizeof(T), DT::MemcpyDeviceToHost);
    for(size_t iB(0);iB<Y.size();++iB)
      ASSERT(B[iB]==iB%10,"Expected content");

    delete[] B;
    return true;
  }

};

typedef Device<CPU> DevCPU;
typedef Device<GPU> DevGPU;

extern const DevCPU cpu_dev(0);
extern const DevGPU gpu_dev(0);

int main(int argc, char* argv[])
{
  VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
  COUT("==============================\n"
       <<" DataIO Test:"
       <<"\n==============================");

  DataIOTest<Array<float, DevCPU, cpu_dev> > io_a_cpu;
  io_a_cpu.PerformAll();

  #ifdef GPU_ACTIVE
      DataIOTest<Array<float, DevGPU, gpu_dev> > io_gpu;
      io_gpu.PerformAll();
  #endif

  VES3D_FINALIZE();
  return 0;
}
