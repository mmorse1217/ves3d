#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include <typeinfo>
#include "enums.h"

template<typename T>
class VesInteraction
{
  public:
    typedef void(*InteractionFun_t)(const T*, const T*, size_t, T*, void*);    

    explicit VesInteraction(InteractionFun_t interaction_handle = NULL, 
        void* user_ptr  = NULL, int num_threads = 1);
    ~VesInteraction();

    template<typename VecContainer>
    InteractionReturn operator()(const VecContainer &position, 
        VecContainer &density, VecContainer &potential) const;

  private:
    InteractionFun_t interaction_handle_;
    int num_threads_;
    void * user_ptr_;    
    size_t* each_thread_np_;
    size_t* each_thread_idx_;

    mutable size_t np_;
    mutable size_t containers_capacity_;

    mutable T* all_pos_;
    mutable T* all_den_;
    mutable T* all_pot_;
    
    size_t getCpyDestIdx(size_t this_thread_np) const;
    void checkContainersSize() const;
    void updatePotential() const;
    
    VesInteraction(VesInteraction const &rhs);
    VesInteraction& operator=(const VesInteraction &rhs);
};

#include "VesInteraction.cc"

#endif // _VESINTERACTION_H_
