#include "CudaApiGlobals.h"

int CudaApiGlobals::num_streams_(0);
cudaStream_t* CudaApiGlobals::streams_(NULL);
int CudaApiGlobals::this_stream_(0);

int CudaApiGlobals::NumStreams()
{
    return(CudaApiGlobals::num_streams_);
}

void CudaApiGlobals::NumStreams(int n)
{
    if ( n == CudaApiGlobals::num_streams_ )
        return;

    COUTDEBUG(" ------------------------------------\n"
        <<"   Creating "<<n<<" CUDA Streams.\n"
        <<" ------------------------------------"<<endl);

    CudaApiGlobals::ClearAll();

    CudaApiGlobals::streams_ = new cudaStream_t[n];
    CudaApiGlobals::num_streams_ = n;
    CudaApiGlobals::this_stream_ = 0;
    
    for (int ii = 0; ii < CudaApiGlobals::num_streams_; ++ii)
        cudaStreamCreate(&CudaApiGlobals::streams_[ii]);
}

const cudaStream_t& CudaApiGlobals::ThisStream()
{
    if( CudaApiGlobals::num_streams_ == 0 )
        CudaApiGlobals::NumStreams(1);

    COUTDEBUG(" ------------------------------------\n"
        <<"   Current Stream "<<CudaApiGlobals::this_stream_<<".\n"
        <<" ------------------------------------"<<endl);

    
    return(CudaApiGlobals::streams_[CudaApiGlobals::this_stream_]);
}

const cudaStream_t& CudaApiGlobals::NextStream()
{
    if( CudaApiGlobals::num_streams_ == 0 )
        CudaApiGlobals::NumStreams(1);
    
    CudaApiGlobals::this_stream_++;
    CudaApiGlobals::this_stream_ %= CudaApiGlobals::num_streams_;

    return( CudaApiGlobals::ThisStream() );
}

void CudaApiGlobals::SyncStream()
{
    for(int ii=0; ii< CudaApiGlobals::num_streams_; ++ii)
        cudaStreamSynchronize(CudaApiGlobals::streams_[ii]);
}

void CudaApiGlobals::SyncStream(int idx)
{
    if (idx < CudaApiGlobals::num_streams_)
        cudaStreamSynchronize(CudaApiGlobals::streams_[idx]);
}

void CudaApiGlobals::ClearAll()
{
    //Destroying all active streams
    for (int ii = 0; ii < CudaApiGlobals::num_streams_; ++ii)
        cudaStreamDestroy(CudaApiGlobals::streams_[ii]);
    
    delete[] CudaApiGlobals::streams_;
    CudaApiGlobals::streams_ = NULL;
    CudaApiGlobals::num_streams_ = 0;
    CudaApiGlobals::this_stream_ = 0;
}
