template<typename T>
Parameters<T> Parameters<T>::sInstance;

template<typename T>
Parameters<T>::Parameters() :
    n_surfs(1),   
    sh_order(12),
    filter_freq(8),
    bending_modulus(1e-2),
    outer_solver_maxit(1),
    inner_solver_maxit(20),
    outer_solver_tol((typeid(T) == typeid(float)) ? 1e-5 : 1e-7),
    inner_solver_tol((typeid(T) == typeid(float)) ? 1e-7 : 1e-9),
    n_steps(1),
    time_horizon(1),
    ts(1),
    rep_maxit(10),
    rep_up_freq(24),
    rep_filter_freq(4),
    rep_ts(1),
    rep_tol(1e-4),
    bg_flow_param(1e-1),
    save_data(false),
    save_stride(-1)
{}

template<typename T>
Parameters<T>::~Parameters()
{}

template<typename T>
const Parameters<T>& Parameters<T>::getInstance()
{
    return sInstance;
}

template<typename T>
Parameters<T>& Parameters<T>::getInstanceModifiable()
{
    return sInstance;
}

template<typename T>
ostream& operator<<(ostream& output, const Parameters<T>& par)
{
    
    output<<"\n ===================================="<<endl;
    output<<"  Simulator parameters"<<endl;
    output<<" ===================================="<<endl;
    output<<"  Surface:"<<endl;
    output<<"    Number of surfaces       : "<<par.n_surfs<<endl;
    output<<"    SH order                 : "<<par.sh_order<<endl;
    output<<"    Filter freq              : "<<par.filter_freq<<endl;
    output<<"    Bending modulus          : "<<par.bending_modulus<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  Solver:"<<endl;
    output<<"    Outer solver maxit       : "<<par.outer_solver_maxit<<endl;
    output<<"    Inner solver maxit       : "<<par.inner_solver_maxit<<endl;
    output<<"    Outer solver tol         : "<<par.outer_solver_tol<<endl;
    output<<"    Inner solver tol         : "<<par.inner_solver_tol<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  Time stepper:"<<endl;
    output<<"    Number of time steps     : "<<par.n_steps<<endl;
    output<<"    Time horizon             : "<<par.time_horizon<<endl;
    output<<"    Step size                : "<<par.ts<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  Reparametrization:"<<endl;
    output<<"    Rep maxit                : "<<par.rep_maxit<<endl;
    output<<"    Rep upsample freq        : "<<par.rep_up_freq<<endl;
    output<<"    Rep filter freq          : "<<par.rep_filter_freq<<endl;
    output<<"    Rep step size            : "<<par.rep_ts<<endl;
    output<<"    Rep tol                  : "<<par.rep_tol<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  Background flow:"<<endl;
    output<<"    Background flow parameter: "<<par.bg_flow_param<<endl;
    output<<" ===================================="<<endl<<endl;
    
    return output;
}
