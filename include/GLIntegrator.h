#ifndef _GLINTEGRATOR_H_
#define _GLINTEGRATOR_H_

#include "HelperFuns.h"

template<typename Container>
class GaussLegendreIntegrator
{
  public:
    GaussLegendreIntegrator();
    ~GaussLegendreIntegrator();

    inline void operator()(const Container &x_in, 
        const Container &w_in, Container &x_dw) const;
    
    inline void operator()(const Container &w_in, Container &dw) const;

  private:
    inline Container* getQuadWeight(int key) const;
    
    DataIO<typename Container::value_type, CPU> IO;
    
    static map<int, Container*> qw_map;
};

#include "GLIntegrator.cc"

#endif _GLINTEGRATOR_H_
