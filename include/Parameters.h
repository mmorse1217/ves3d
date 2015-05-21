#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include "Enums.h"
#include "Error.h"
#include "Logger.h"
#include "anyoption.h"
#include "Streamable.h"
#include "ves3d_common.h"

#include <typeinfo>
#include <iostream>
#include <string>
#include <map>


typedef std::map<std::string,std::string> DictString_t;

template <typename T>
class Parameters : public Streamable
{
  public:
    Parameters();
    Parameters(int argc, char** argv);
    Parameters(std::istream &is, Format format);

    ~Parameters();

    //Surface
    int n_surfs;
    int sh_order;
    int filter_freq;
    int upsample_freq;

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
    bool time_adaptive;
    bool solve_for_velocity;
    bool pseudospectral;

    enum SolverScheme scheme;
    enum PrecondScheme time_precond;
    enum BgFlowType bg_flow;
    enum SingularStokesRot singular_stokes;

    //Reparametrization
    int  rep_maxit;
    int  rep_filter_freq;
    bool rep_upsample;
    T    rep_ts;
    T    rep_tol;

    //Background flow
    T bg_flow_param;
    bool interaction_upsample;
    T periodic_length;

    //Startup and Monitoring
    bool checkpoint;
    T checkpoint_stride;
    std::string init_file_name;
    std::string cntrs_file_name;
    std::string checkpoint_file_name;
    std::string load_checkpoint;
    T error_factor;
    int num_threads;

    //parsing
    Error_t parseInput(int argc, char** argv, const DictString_t *dict=NULL);
    Error_t expand_templates(const DictString_t *dict=NULL);

    //Adjust the frequencies according to the sh_order
    void adjustFreqs();

    // From streamable class --------------------------------------------------
    // ------------------------------------------------------------------------
    virtual Error_t pack(std::ostream &os, Format format) const;
    virtual Error_t unpack(std::istream &is, Format format);

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

//! Expands the templates in pattern marked by their key in dict to
//! the value in dict. For example user_{{name}}} with dict
//! {name:'JohnDoe'} expands to user_JohnDoe.
Error_t expand_template(std::string *pattern, const DictString_t &dict);

#include "Parameters.cc"

#endif //_PARAMETERS_H_
