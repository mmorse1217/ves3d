#ifndef _CUDAAPIGLOBALS_H_
#define _CUDAAPIGLOBALS_H_

#include "cuda_runtime.h"
#include "Logger.h"

class CudaApiGlobals
{
  public:  
    static int NumStreams();
    static void NumStreams(int n);

    static const cudaStream_t& ThisStream();
    static const cudaStream_t& NextStream();

    static void SyncStream();
    static void SyncStream(int idx);

    static void ClearAll();
    
  private:
    CudaApiGlobals();
    static int num_streams_;
    static cudaStream_t* streams_;
    static int this_stream_;
};

#endif //_CUDAAPIGLOBALS_H_
