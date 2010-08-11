#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

#include <omp.h>

/**
 * The abstract base class for the repartitioning. Derived classes
 * will act as a gateway to any external repartitioning
 * library/function available. The typedef <tt> GlobalRepart_t</tt> is
 * the expected signature of the external function.
 */
template<typename T>
class RepartitionBase
{
  public:
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride, 
        const T* x, const T* tension, size_t* nvr, T** xr,
        T** tensionr, void* user_ptr);
    
    ~RepartitionBase() = 0;

    template<typename Container>
    void operator()(Container &coord, Container &tension, 
        void* user_ptr = NULL) const;
  
  private:
    RepartitionBase(RepartitionBase const &rhs);
    RepartitionBase& operator=(const RepartitionBase &rhs);
};


template<typename T>
class Repartition : public RepartitionBase<T>
{
  public:
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride, 
        const T* x, const T* tension, size_t* nvr, T** xr,
        T** tensionr, void* user_ptr);
    
    explicit Repartition(GlobalRepart_t fun_ptr = NULL, 
        int num_threads = omp_get_max_threads());
    ~Repartition();

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
};

#include "Repartition.cc"

#endif // _REPARTIONGATEWAY_H_
