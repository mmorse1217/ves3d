#ifndef _EVALVELOCITY_H_
#define _EVALVELOCITY_H_

tempalate<typename Container>
class EvalVelocity
{
  public:
    EvalVelocity();
    ~EvalVelocity();

    Error_t operator()(const Container &x_src, 
        const Container &x_eval, Container &vel);

  private:

};

#include "EvalVelocity.cc"

#endif //_EVALVELOCITY_H_
