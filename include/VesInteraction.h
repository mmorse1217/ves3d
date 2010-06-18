#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include "HelperFuns.h"

template<typename VecContainer>
class VesInteraction
{
  public:
    void operator()(const VecContainer &position, VecContainer &density, 
        VecContainer &potential) const 
    {
        VecContainer wrk1, wrk2, wrk3;
        wrk1.replicate(position);
        wrk2.replicate(position);
        wrk3.replicate(position);
        
        size_t np(position.getNumSubs() * position.getStride());
        
        position.getDevice().Transpose(position.begin(), np, 3, wrk1.begin());
        position.getDevice().Transpose( density.begin(), np, 3, wrk2.begin());
        
        position.getDevice().DirectStokes(wrk1.begin(), wrk2.begin(), NULL, 
            np, 1, wrk1.begin(), 0, np, wrk3.begin());  
        
        position.getDevice().Transpose(wrk3.begin(), 3, np, potential.begin());
    }
};

#endif // _VESINTERACTION_H_
