#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

#include <omp.h>
#include <cassert>
#include "Error.h"

/**
 * The interface class with the global repartitioning function. It
 * gathers data from multiple threads and passed them to the external
 * repartitioning function, which repartitions data to match the MPI
 * load balance. After returning form the external function, this
 * class distributes the vesicles between the available threads.
 */
template<typename T>
class Repartition
{
  public:
    ///The function pointer type for the external repartitioning function.
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride,
        const T* x, const T* tension, size_t* nvr, T** xr,
        T** tensionr, void** context);

    //Deallocator for the context
    typedef void(*Dealloc_t)(void**);

    /**
     * @param fun_ptr The actual pointer to the repartitioning
     * function, no repartitioning when <tt>NULL</tt>.
     * @param num_threads The number of expected threads that to be handled.
     */
    explicit Repartition(GlobalRepart_t fun_ptr = NULL,
        Dealloc_t clear_context = NULL,
        int num_threads = omp_get_max_threads());
    ~Repartition();

    /**
     * The function called from each thread.
     * @param coord The Cartesian coordinate of the points
     * @param tension The tension associated with each point.
     * @param user_ptr the user-defined pointer that may be needed
     * depending on the external repartitioning function.
     */
    template<typename VecContainer, typename ScaContainer>
    Error_t operator()(VecContainer &coord, ScaContainer &tension) const;

  private:
    GlobalRepart_t g_repart_handle_;
    Dealloc_t clear_context_;

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
    mutable void* context_;

    size_t getCpyIdx(size_t this_thread_nv, size_t stride) const;
    size_t getNvShare() const;
    void checkContainersSize(size_t stride) const;
};

#include "Repartition.cc"

#endif // _REPARTIONGATEWAY_H_
