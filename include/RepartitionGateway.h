#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

#include <omp.h>

template<typename Container>
class RepartitionGateway
{
  private:
    typedef typename Container::value_type value_type;
    
  public:
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride, 
        const value_type* x, const value_type* tension, size_t* nvr, 
        value_type** xr, value_type** tensionr, void* user_ptr);
    
    explicit RepartitionGateway(GlobalRepart_t fun_ptr = NULL, 
        int num_threads = 1);
    ~RepartitionGateway();

    void operator()(Container &coord, Container &tension, 
        void* user_ptr) const;
    
  private:
    GlobalRepart_t g_repart_handle_;
    int num_threads_;
    size_t* each_thread_nv_;
    size_t* each_thread_idx_;

    mutable size_t nv_;
    mutable size_t capacity_;
    
    mutable value_type* all_pos_;
    mutable value_type* all_tension_;
    
    mutable value_type* posr_;
    mutable value_type* tensionr_;
    mutable size_t nvr_;


    size_t getCpyIdx(size_t this_thread_nv, size_t stride) const;
    size_t getNvShare() const;
    void checkContainersSize() const;
        
    RepartitionGateway(RepartitionGateway const &rhs);
    RepartitionGateway& operator=(const RepartitionGateway &rhs);
};

#include "RepartitionGateway.cc"

#endif // _REPARTIONGATEWAY_H_
