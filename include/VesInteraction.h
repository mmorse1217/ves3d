#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include<typeinfo>

template<typename T>
class VesInteraction
{
  public:
    enum InteractionReturn {InteractionSuccess = 0,
                            NoInteraction      = 1};

    typedef void(*InteractionFun)(const T*, const T*, size_t, T*);    

    VesInteraction(InteractionFun interaction_handle = NULL, int num_threads = 1);
    ~VesInteraction();

    template<typename VecContainer>
    InteractionReturn operator()(const VecContainer &position, VecContainer &density,
        VecContainer &potential) const;

  private:
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
    void updatePotential() const;
    InteractionFun interaction_handle_;
    
    VesInteraction(VesInteraction const &rhs);
    VesInteraction& operator=(const VesInteraction &rhs);
};

#include "VesInteraction.cc"

#endif // _VESINTERACTION_H_
