#ifndef _GLINTEGRATOR_H_
#define _GLINTEGRATOR_H_

#include "HelperFuns.h"
#include "VesBlas.h"
#include <typeinfo>
#include "DataIO.h"
#include "Spharm.h"

template<typename Container>
class GaussLegendreIntegrator
{
  public:
    GaussLegendreIntegrator();
    ~GaussLegendreIntegrator();

    template<typename InputContainer>
    inline void operator()(const InputContainer &x_in,
        const Container &w_in, InputContainer &x_dw) const;

    inline void operator()(const Container &w_in,
        Container &dw) const;

    Container* getQuadWeights(int key) const;
    Container* buildQuadWeights(int shOrder) const;

  private:
    static std::map<int, Container*> qw_map;
};

#include "GLIntegrator.cc"

#endif //_GLINTEGRATOR_H_
