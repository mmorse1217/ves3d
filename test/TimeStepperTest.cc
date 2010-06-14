#include "TimeStepper.h"
#include <iostream>
#include <math.h>
using namespace std;

#ifndef Doxygen_skip

typedef float real;

class CC
{
  public:
    typedef real value_type;
    CC(real val_in = 0) : val(val_in) {}; 
    real val;
};

class FC
{
  public:
    virtual void  operator()(const CC &state, 
        const real &t, CC &val) const 
    {
        val.val = -2*state.val;
    }   
};

class MC 
{
  public:
    MC(real hor){ time_horizon_ = hor;}
    
    virtual bool operator()(const CC &state, 
        const real &t, real &dt) const
    {
        cout<<state.val-exp(-2*t)<<endl;
        return(t<time_horizon_);
    }
    
    real time_horizon_;
};

class ForwardEuler
{
  public:
    virtual void operator()(const CC &old_state, 
        const real &t, const real &dt, 
        const FC &F, CC &new_state) const
    {
        F(old_state.val, t, temp);
        new_state.val = old_state.val + dt * temp.val;
    }

    mutable CC temp;
};

#endif //Doxygen_skip

int main(int argc, char **argv)
{

    CC y(1);
    CC::value_type t = 0;
    CC::value_type dt = .001;
    FC Fun;
    ForwardEuler Dis;
    MC Mntr(1);
    
    TimeStepper<CC,FC, ForwardEuler, MC> Stepper;
    Stepper(y, t, dt, Fun, Dis, Mntr);

    return 0;
}
