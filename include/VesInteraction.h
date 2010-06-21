#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

//#include "HelperFuns.h"

template<typename VecContainer>
class VesInteraction
{
  public:
    typedef typename VecContainer::value_type value_type;
    typedef double fmm_value_type;
    
    VesInteraction(int num_threads);
    ~VesInteraction();

    void operator()(const VecContainer &position, VecContainer &density,
        VecContainer &potential);

  private:
    int num_threads_;
    size_t* each_thread_np_;
    size_t* each_thread_idx_;
    
    size_t np_;
    size_t containers_capacity_;

    fmm_value_type* all_pos_;
    fmm_value_type* all_den_;
    fmm_value_type* all_pot_;
    
    size_t getCpyDestIdx(size_t this_thread_np);
    void checkContainersSize();

    void fmmInteraction() const;
    void directInteraction() const;

    VesInteraction(VesInteraction const &rhs);
    VesInteraction& operator=(const VesInteraction &rhs);
};

#include "VesInteraction.cc"

#endif // _VESINTERACTION_H_
