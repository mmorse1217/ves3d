template<typename T>
Parameters<T>::Parameters() :
    n_surfs(1),   
    sh_order(12),
    filter_freq(8),
    bending_modulus(1e-2),
    position_solver_iter(15),
    tension_solver_iter(15),
    position_solver_restart(1),
    tension_solver_restart(1),
    position_solver_tol((typeid(T) == typeid(float)) ? 1e-4: 1e-8),
    tension_solver_tol((typeid(T) == typeid(float)) ? 5e-4: 1e-8),
    n_steps(1),
    time_horizon(1),
    ts(1),
    scheme(Explicit),
    singular_stokes(ViaSpHarm),
    rep_maxit(10),
    rep_up_freq(24),
    rep_filter_freq(4),
    rep_ts(1),
    rep_tol((typeid(T) == typeid(float)) ? 1e-3 : 1e-4),
    bg_flow_param(1e-1),
    save_data(false),
    save_stride(-1),
    error_factor(1)
{}

template<typename T>
Parameters<T>::~Parameters()
{}

template<typename T>
std::ostream& operator<<(std::ostream& output, const Parameters<T>& par)
{
    output<<"\n ===================================="<<std::endl;
    output<<"  Simulator parameters"<<std::endl;
    output<<" ===================================="<<std::endl;
    output<<"  Surface:"<<std::endl;
    output<<"    Number of surfaces       : "<<par.n_surfs<<std::endl;
    output<<"    SH order                 : "<<par.sh_order<<std::endl;
    output<<"    Filter freq              : "<<par.filter_freq<<std::endl;
    output<<"    Bending modulus          : "<<par.bending_modulus<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Solver:"<<std::endl;
    output<<"    Position solver iter     : "<<par.position_solver_iter<<std::endl;
    output<<"    tension solver iter      : "<<par.tension_solver_iter<<std::endl;
    output<<"    Position solver restart  : "<<par.position_solver_restart<<std::endl;
    output<<"    Tension solver restart   : "<<par.position_solver_restart<<std::endl;
    output<<"    Position solver tol      : "<<par.position_solver_tol<<std::endl;
    output<<"    Tension solver tol       : "<<par.tension_solver_tol<<std::endl;

    output<<" ------------------------------------"<<std::endl;
    output<<"  Time stepper:"<<std::endl;
    output<<"    Number of time steps     : "<<par.n_steps<<std::endl;
    output<<"    Time horizon             : "<<par.time_horizon<<std::endl;
    output<<"    Step size                : "<<par.ts<<std::endl;
    output<<"    Scheme                   : "<<par.scheme<<std::endl;
    output<<"    Singular Stokes          : "<<par.singular_stokes<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Reparametrization:"<<std::endl;
    output<<"    Rep maxit                : "<<par.rep_maxit<<std::endl;
    output<<"    Rep upsample freq        : "<<par.rep_up_freq<<std::endl;
    output<<"    Rep filter freq          : "<<par.rep_filter_freq<<std::endl;
    output<<"    Rep step size            : "<<par.rep_ts<<std::endl;
    output<<"    Rep tol                  : "<<par.rep_tol<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Saving:"<<std::endl;
    output<<"    Save data                : "<<std::boolalpha<<par.save_data<<std::endl;
    output<<"    Save file name           : "<<par.save_file_name<<std::endl;
    output<<"    Save stride              : "<<par.save_stride<<std::endl;
    output<<"    Error Factor             : "<<par.error_factor<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Background flow:"<<std::endl;
    output<<"    Background flow parameter: "<<par.bg_flow_param<<std::endl;
    output<<" ===================================="<<std::endl<<std::endl;

    return output;
}
