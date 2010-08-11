#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include <typeinfo>
#include "enums.h"

/**
 * The abstract base class for interaction between vesicles. Derived
 * classed act as a gateway between the local code and the far
 * interaction FMM code. They take care of multi-threads, orgainzing
 * date, copying to the host, etc. The typedef
 * <tt>ExternalInteractionFun_t</tt> defines the expected signature of
 * the FMM function.
 */
template<typename T>
class VesInteractionBase
{
  public:
    typedef void(*InteractionFun_t)(const T* pos, const T*den, size_t np, 
        T* pot, void* usr_ptr);    
    
    ~VesInteractionBase() = 0;

    template<typename VecContainer>
    InteractionReturn operator()(const VecContainer &position, 
        VecContainer &density, VecContainer &potential,
        void* usr_ptr = NULL) const;
};

///The gateway function between the local code and FMM.
template<typename T>
class VesInteraction : public VesInteractionBase<T>
{
  public:
    typedef void(*InteractionFun_t)(const T*, const T*, size_t, T*, void*);    
    
    explicit VesInteraction(InteractionFun_t interaction_handle = NULL,
        int num_threads = omp_get_max_threads());
    ~VesInteraction();

    template<typename VecContainer>
    InteractionReturn operator()(const VecContainer &position, 
        VecContainer &density, VecContainer &potential,
        void* usr_ptr = NULL) const;

  private:
    InteractionFun_t interaction_handle_;
    int num_threads_;
    size_t* each_thread_np_;
    size_t* each_thread_idx_;

    mutable size_t np_;
    mutable size_t containers_capacity_;

    mutable T* all_pos_;
    mutable T* all_den_;
    mutable T* all_pot_;
    
    size_t getCpyDestIdx(size_t this_thread_np) const;
    void checkContainersSize() const;
    
    VesInteraction(VesInteraction const &rhs);
    VesInteraction& operator=(const VesInteraction &rhs);
};

#include "VesInteraction.cc"

#endif // _VESINTERACTION_H_
