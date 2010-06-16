#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include<typeinfo>

template <typename T>
struct Parameters 
{
  public:
    //Surface
    int n_surfs;   
    int sh_order;
    int filter_freq;
    T bending_modulus;
    
    //Solver parameters
    int outer_solver_maxit;
    int inner_solver_maxit;
    T outer_solver_tol;
    T inner_solver_tol;

    //Time stepper
    int n_steps;
    T time_horizon;
    T ts;
    
    //Reparametrization
    int rep_maxit;
    int rep_up_freq;
    int rep_filter_freq;
    T rep_ts;
    T rep_tol;

    //Background flow
    T bg_flow_param;

    //Monitoring
    bool save_data;
    int save_stride;
    
    //Accessor method   
    static const Parameters& getInstance();
    static Parameters& getInstanceModifiable();

  protected:
    Parameters();
    ~Parameters();
    
    static Parameters<T> sInstance;
};

template<typename T>
ostream& operator<<(ostream& output, const Parameters<T>& par);

#include "Parameters.cc"

#endif //_PARAMETERS_H_
