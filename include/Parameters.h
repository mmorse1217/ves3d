#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <typeinfo>
#include <iostream>
#include <string>
#include "Enums.h"
#include "Error.h"
#include "Logger.h"
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
  T viscosity_contrast;

  //Solver parameters
  int position_solver_iter;
  int tension_solver_iter;
  int position_solver_restart;
  int tension_solver_restart;
  T position_solver_tol;
  T tension_solver_tol;

  //Time stepper
  T time_horizon;
  T ts;
  T time_tol;
  int time_iter_max;

  enum SolverScheme scheme;
  enum PrecondScheme time_precond;
  enum SingularStokesRot singular_stokes;

  //Reparametrization
  int rep_maxit;
  int rep_up_freq;
  int rep_filter_freq;
  T rep_ts;
  T rep_tol;

  //Background flow
  T bg_flow_param;
  bool upsample_interaction;

  //Startup and Monitoring
  bool save_data;
  T save_stride;
  std::string init_file_name;
  std::string cntrs_file_name;
  std::string save_file_name;
  T error_factor;
  int num_threads;

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
  void init();
};

template<typename T>
std::ostream& operator<<(std::ostream& output, const Parameters<T>& par);

#include "Parameters.cc"

#endif //_PARAMETERS_H_
