#include "TimeStepper.h"
#include <iostream>

using namespace std;

class BB
{
  public:
    typedef float value_type;
    BB(float val_in = 0){ val = val_in;}
    float val;
};

class CC
{
  public:
    typedef float value_type;
    CC(float val_in = 0) : val(val_in) {}; 
    BB val;
};

class FC : public Forcing<BB>
{
  public:
    virtual void  operator()(const BB &state, 
        const float &t, BB &val) const 
    {
        val.val = -state.val;
    }   
};

class MC : public Monitor<CC>
{
  public:
    MC(float hor){ time_horizon_ = hor;}
    
    virtual bool operator()(const CC &state, 
        const float &t, float &dt) const
    {
        cout<<state.val.val<<endl;
        return(t<time_horizon_);
    }
    
    float time_horizon_;
};

class ForwardEuler : public Discretization<CC, FC>
{
  public:
    virtual void operator()(const CC &old_state, 
        const float &t, const FC &F, 
        const float &dt, CC &new_state) const
    {
        F(old_state.val, t, temp);
        new_state.val.val = old_state.val.val + dt * temp.val;
    }

    mutable BB temp;
};

int main(int argc, char **argv)
{

    CC y(1);
    CC::value_type t = 0;
    CC::value_type dt = .0001;
    FC Fun;
    ForwardEuler Dis;
    MC Mntr(1);
    
    TimeStepper<CC,FC> Stepper;
    Stepper(y, t, dt, Fun, Dis, Mntr);

    return 0;
}
