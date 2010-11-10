#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <typeinfo>
#include <iostream>
#include <string>
#include "enums.h"
#include "anyoption.h"

template <typename T>
struct Parameters 
{
public:
  Parameters();
  Parameters(int argc, char** argv);

  ~Parameters();

  //Surface
  int n_surfs;   
  int sh_order;
  int filter_freq;
  T bending_modulus;
    
  //Solver parameters
  int position_solver_iter;
  int tension_solver_iter;
  int position_solver_restart;
  int tension_solver_restart;
  T position_solver_tol;
  T tension_solver_tol;

  //Time stepper
  int n_steps;
  T time_horizon;
  T ts;
  enum SolverScheme scheme;
  enum SingularStokesRot singular_stokes;

  //Reparametrization
  int rep_maxit;
  int rep_up_freq;
  int rep_filter_freq;
  T rep_ts;
  T rep_tol;  

  //Background flow
  T bg_flow_param;

  //Startup and Monitoring
  bool save_data;
  T save_stride;
  std::string init_file_name;
  std::string cntrs_file_name;
  std::string save_file_name;
  T error_factor;

  //parsing
  Error_t parseInput(int argc, char** argv);

  //Adjust the frequencies according to the sh_order
  void adjustFreqs();

private:
  Parameters(Parameters<T> &rhs);
  Parameters<T>& operator=(Parameters<T> &rhs);

  void setUsage(AnyOption* opt);
  void setOptions(AnyOption* opt);
  void getOptionValues(AnyOption* opt);
};

template<typename T>
std::ostream& operator<<(std::ostream& output, const Parameters<T>& par);

#include "Parameters.cc"

#endif //_PARAMETERS_H_
