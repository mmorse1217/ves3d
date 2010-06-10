#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

template<typename Container,
         typename Forcing,
         typename Updater,
         typename Monitor>
class TimeStepper
{
  private:
    typedef typename Container::value_type value_type;

  public:
    void operator()(Container &state, value_type &t, value_type &dt,
        Forcing &F, Updater &U, Monitor &M)
    {
        while ( M(state, t, dt) )
        {
            U(state, t, dt, F, state);
            t += dt;
        }
    }
};

#endif //_TIMESTEPPER_H_
