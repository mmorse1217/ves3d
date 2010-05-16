// The parameter struct methods
template<typename T>
SurfaceParams<T>::SurfaceParams()
{
    this->mapStringValues["p"]               = 1;
    this->mapStringValues["n_surfs"]         = 2;
    this->mapStringValues["kappa"]           = 3;
    this->mapStringValues["filter_freq"]     = 4;
    this->mapStringValues["rep_ts"]          = 5;
    this->mapStringValues["rep_max_vel"]     = 6;
    this->mapStringValues["rep_iter_max"]    = 7;
    this->mapStringValues["rep_up_freq"]     = 8;
    this->mapStringValues["rep_filter_freq"] = 9;
}

template<typename T>
void SurfaceParams<T>::SetMember(string var_name, string var_val)
{
    char *cstr = new char [var_val.size()+1];
    float val_num;
    
    switch (this->mapStringValues[var_name])
    {
        case 1:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->p_ = val_num;
            break;

        case 2:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->n_surfs_ = val_num;
            break;

        case 3:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->kappa_ = val_num;
            break;

        case 4:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->filter_freq_ = val_num;
            break;

        case 5:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->rep_ts_ = val_num;
            break;

        case 6:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->rep_max_vel_ = val_num;
            break;

        case 7:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->rep_iter_max_ = val_num;
            break;
    
        case 8:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->rep_up_freq_ = val_num;
            break;

        case 9:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->rep_filter_freq_ = val_num;
            break;

        default:
            cout << "'" << var_name
                 << "' is an invalid parameter."
                 << endl;
            
            break;
    }
    delete cstr;
}

template<typename T>
ostream& operator<<(ostream& output, const SurfaceParams<T>& par)
{
    
    output<<"\n ------------------------------------"<<endl;
    output<<"  Surface properties"<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  p                 : "<<par.p_<<endl;
    output<<"  Number of surfaces: "<<par.n_surfs_<<endl;
    output<<"  kappa             : "<<par.kappa_<<endl;
    output<<"  filter_freq       : "<<par.filter_freq_<<endl;
    output<<"  rep_ts            : "<<par.rep_ts_<<endl;
    output<<"  rep_max_iter      : "<<par.rep_iter_max_<<endl;
    output<<"  rep_max_vel       : "<<par.rep_max_vel_<<endl;
    output<<"  rep_up_freq       : "<<par.rep_up_freq_<<endl;
    output<<"  rep_filter_freq   : "<<par.rep_filter_freq_<<endl;
    output<<" ------------------------------------"<<endl<<endl;
    
    return output;
}
