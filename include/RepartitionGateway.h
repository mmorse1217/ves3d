#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

#include <omp.h>

template<typename T>
class RepartitionGateway
{
  public:
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride, 
        const T* x, const T* tension, size_t* nvr, T** xr,
        T** tensionr, void* user_ptr);
    
    explicit RepartitionGateway(GlobalRepart_t fun_ptr = NULL, 
        int num_threads = 1);
    ~RepartitionGateway();

    template<typename Container>
    void operator()(Container &coord, Container &tension, 
        void* user_ptr = NULL) const;
    
  private:
    GlobalRepart_t g_repart_handle_;
    int num_threads_;
    size_t* each_thread_nv_;
    size_t* each_thread_idx_;

    mutable size_t nv_;
    mutable size_t capacity_;
    
    mutable T* all_pos_;
    mutable T* all_tension_;
    
    mutable T* posr_;
    mutable T* tensionr_;
    mutable size_t nvr_;


    size_t getCpyIdx(size_t this_thread_nv, size_t stride) const;
    size_t getNvShare() const;
    void checkContainersSize(size_t stride) const;
        
    RepartitionGateway(RepartitionGateway const &rhs);
    RepartitionGateway& operator=(const RepartitionGateway &rhs);
};

#include "RepartitionGateway.cc"

#endif // _REPARTIONGATEWAY_H_
