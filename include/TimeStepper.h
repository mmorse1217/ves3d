#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

template<typename Container>
class Forcing
{
  public:
    virtual void  operator()(const Container &state, 
        const typename Container::value_type &t, 
        Container &val) const = 0;
};

template<typename Container>
class Monitor
{
  public:
    virtual bool operator()(const Container &state, 
        const typename Container::value_type &t, 
        typename Container::value_type &dt) const = 0;
};

template<typename Container, typename Forcing>
class Discretization
{
  public:
    virtual void operator()(const Container &old_state, 
        const typename Container::value_type &t, 
        const Forcing &F, const typename Container::value_type &dt, 
        Container &new_state) const = 0;
};

template<typename Container, typename Forcing>
class TimeStepper
{
  public:
    typedef typename Container::value_type value_type;

    void operator()(Container &state, 
        typename Container::value_type &time,
        typename Container::value_type &dt, Forcing &F, 
        Discretization<Container, Forcing> &D, 
        Monitor<Container> &M)
    {
        while ( M(state, time, dt) )
        {
            D(state, time, F, dt, state);
            time+=dt;
        }
    }
};

#endif //_TIMESTEPPER_H_
