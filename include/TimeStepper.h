#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"
#include "VesUtil.h"
#include "DataIO.h"

#define PI_8I 1.0/8.0/M_PI

template<typename T> 
class TimeStepper
{
  public:
    T ts_;
    int n_steps_;
    DataIO<T> &fileIO;
    Surface<T> &vesicle_;
    Vectors<T> velocity_, vel_bending_, vel_tension_;

    VelField<T> &bg_flow_;
    void(*Interaction_)(T*, T*, int, int, T*);

    T *quad_weights_;

    TimeStepper(T ts_in, int n_steps_in, Surface<T> &ves_in, 
        DataIO<T> &fileIO_in, VelField<T> &bg_flow_in, 
        void(*Interaction_in)(T*, T*, int, int,  T*)) : 
        ts_(ts_in),
        n_steps_(n_steps_in), 
        fileIO(fileIO_in),
        vesicle_(ves_in),
        velocity_   (vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_bending_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        vel_tension_(vesicle_.device_, vesicle_.params_.p_, vesicle_.params_.n_surfs_),
        bg_flow_(bg_flow_in),
        Interaction_(Interaction_in)
    {
        int np = 2 * vesicle_.params_.rep_up_freq_ * (vesicle_.params_.rep_up_freq_ + 1);
        quad_weights_ = (T*) malloc(np * sizeof(T));
        char fname[300];
        sprintf(fname,"%s/precomputed/quad_weights_%u_single.txt",getenv("VES3D_DIR"),vesicle_.params_.rep_up_freq_);
        fileIO.ReadData(fname, np, quad_weights_);
    };

    ~TimeStepper()
    {
        free(quad_weights_);
    };

    void EvolveInTime()
    {
        cout<<vesicle_.params_;
        
        cout<<" ------------------------------------"<<endl;
        cout<<"  Time stepper parameters"<<endl;
        cout<<" ------------------------------------"<<endl;
        cout<<"  ts                : "<<ts_<<endl;
        cout<<"  n_steps           : "<<n_steps_<<endl;
        cout<<" ------------------------------------"<<endl<<endl;

        for(int idx=0;idx<n_steps_;idx++)
        {

            cout<<idx<<endl;
            vesicle_.UpdateAll();

            //Interaction
            avpw(vesicle_.tension_, vesicle_.tensile_force_, vesicle_.bending_force_, vel_tension_);
            //vel_tension_ holds the interaction force
            xvpb(vesicle_.w_,vel_tension_, (T) 0.0, vel_tension_);
            axpb((T) (PI_8I), vel_tension_, (T) 0.0, vel_tension_);

            //up-sampling x_ to V10
            vesicle_.device_.ShAna(vesicle_.x_.data_, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V10.data_);
            vesicle_.device_.Resample(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, vesicle_.params_.rep_up_freq_, 
                vesicle_.V10.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V10.data_);

            //up-sampling vel_tension to V11
            vesicle_.device_.ShAna(vel_tension_.data_, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);
            vesicle_.device_.Resample(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, vesicle_.params_.rep_up_freq_, 
                vesicle_.V11.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);

            //the self-interaction -- V12
            vesicle_.device_.DirectStokes(vesicle_.V10.GetFunLength(), vesicle_.params_.n_surfs_, 
                0, vesicle_.V10.GetFunLength(), quad_weights_, vesicle_.V10.data_,
                vesicle_.V10.data_, vesicle_.V11.data_, vesicle_.V12.data_);

            ///@todo the multiplication by the quadrature weights is slow
            int stride = vesicle_.V11.GetFunLength();
            for(int ii=0;ii<vesicle_.params_.n_surfs_;++ii)
            {
                vesicle_.device_.xvpb(quad_weights_, vesicle_.V11.data_ + 3*ii*stride,
                    (T) 0.0, stride, 1, vesicle_.V11.data_ + 3*ii*stride);
            }

            //Interaction to V13
            Interaction_(vesicle_.V10.data_, vesicle_.V11.data_,  
                vesicle_.V10.GetFunLength(), vesicle_.params_.n_surfs_, vesicle_.V13.data_);
            
            axpy((T) -1.0, vesicle_.V12, vesicle_.V13, vesicle_.V13);
            
            vesicle_.device_.ShAna(vesicle_.V13.data_, vesicle_.work_arr, vesicle_.params_.rep_up_freq_, 
                3*vesicle_.params_.n_surfs_, vesicle_.V11.data_);
            
            //filtering to 1/3
            vesicle_.device_.ScaleFreqs(vesicle_.params_.rep_up_freq_, 3*vesicle_.params_.n_surfs_, 
                vesicle_.V11.data_, vesicle_.alpha_q, vesicle_.V11.data_);

            vesicle_.device_.Resample(vesicle_.params_.rep_up_freq_, 3*vesicle_.params_.n_surfs_, 
                vesicle_.params_.p_, vesicle_.V11.data_, vesicle_.shc);
            vesicle_.device_.ShSyn(vesicle_.shc, vesicle_.work_arr, vesicle_.params_.p_, 
                3*vesicle_.params_.n_surfs_, velocity_.data_);

            //Filtering
            //vesicle_.device_.Filter(vesicle_.params_.p_, 3*vesicle_.params_.n_surfs_, 
            //    velocity_.data_, vesicle_.alpha_p, vesicle_.work_arr, vesicle_.shc, velocity_.data_);

            //DotProduct(velocity_, velocity_, vesicle_.S1);
            //cout<<vesicle_.S1.Max()<<endl;
            //Background flow
            bg_flow_.GetVel(vesicle_.x_, vel_tension_);
            axpy((T) 1.0, vel_tension_, velocity_, velocity_);
            
            //Calculate stokes
            vesicle_.StokesMatVec(vesicle_.bending_force_, vel_bending_);
            axpy((T) 1.0, vel_bending_, velocity_, velocity_);

            vesicle_.StokesMatVec(vesicle_.tensile_force_, vel_tension_);

            //Calculate tension
            vesicle_.GetTension(velocity_, vel_tension_, vesicle_.tension_);
            vesicle_.device_.avpw(vesicle_.tension_, vel_tension_.data_,
                velocity_.data_, velocity_.GetFunLength(), vesicle_.params_.n_surfs_, velocity_.data_);

            //Advance in time
            axpy(ts_, velocity_, vesicle_.x_, vesicle_.x_);
        
            vesicle_.Reparam();

            if((100*(idx+1))%n_steps_ == 0)
                 fileIO.Append(vesicle_.x_.data_, vesicle_.x_.GetDataLength());
        }
    };
};


///All parameters class ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct TimeStepperParams 
{
    T ts_;
    int n_steps_;
    T shear_rate_;

    TimeStepperParams();
    void SetMember(string var_name, string var_val);
    
  private:
    map<string, int> mapStringValues;
};

template<typename T>
TimeStepperParams<T>::TimeStepperParams()
{
    this->mapStringValues["ts"]         = 1;
    this->mapStringValues["n_steps"]    = 2;
    this->mapStringValues["shear_rate"] = 3;
}

template<typename T>
void TimeStepperParams<T>::SetMember(string var_name, string var_val)
{
    char *cstr = new char [var_val.size()+1];
    float val_num;
    
    switch (this->mapStringValues[var_name])
    {
        case 1:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->ts_ = val_num;
            break;

        case 2:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->n_steps_ = val_num;
            break;

        case 3:
            strcpy(cstr, var_val.c_str());
            sscanf(cstr," %g",&val_num);
            this->shear_rate_ = val_num;
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
ostream& operator<<(ostream& output, const TimeStepperParams<T>& par)
{
    
    output<<"\n ------------------------------------"<<endl;
    output<<"  Time Stepper properties"<<endl;
    output<<" ------------------------------------"<<endl;
    output<<"  ts              : "<<par.ts_<<endl;
    output<<"  shear_rate      : "<<par.shear_rate_<<endl;
    output<<"  Number of steps : "<<par.n_steps_<<endl;
    output<<" ------------------------------------"<<endl<<endl;
    
    return output;
}    

struct FileList
{
    string initial_shape_file_;        
    string simulation_out_file_;
    string leg_transform_file_;        
    string inv_leg_transform_file_;    
    string inv_d1_leg_transform_file_; 
    string inv_d2_leg_transform_file_; 
    string all_rotation_matrices_file_;
    string signualar_quad_weight_file_;
    string regular_quad_weigh_file_;   
    string w_sph_file_;                 
    
    FileList();
    void SetMember(string var_name, string var_val);
    
  private:
    map<string, int> mapStringValues;

};

FileList::FileList()
{
    this->mapStringValues["initial_shape_file"]         = 1;        
    this->mapStringValues["simulation_out_file"]        = 2;
    this->mapStringValues["leg_transform_file"]         = 3;        
    this->mapStringValues["inv_leg_transform_file"]     = 4;    
    this->mapStringValues["inv_d1_leg_transform_file"]  = 5; 
    this->mapStringValues["inv_d2_leg_transform_file"]  = 6; 
    this->mapStringValues["all_rotation_matrices_file"] = 7;
    this->mapStringValues["signualar_quad_weight_file"] = 8;
    this->mapStringValues["regular_quad_weigh_file"]    = 9;   
    this->mapStringValues["w_sph_file"]                 = 10;                

}

void FileList::SetMember(string var_name, string var_val)
{    
    switch (this->mapStringValues[var_name])
    {
        case 1:
            this->initial_shape_file_ = var_val;   
            break;

        case 2:
            simulation_out_file_ = var_val;
            break;

        case 3:
            leg_transform_file_ = var_val; 
            break;        

        case 4:
            inv_leg_transform_file_ = var_val;
            break;    
            
        case 5:
            inv_d1_leg_transform_file_ = var_val;
            break;

        case 6:
            inv_d2_leg_transform_file_ = var_val;
            break; 
            
        case 7:
            all_rotation_matrices_file_ = var_val; 
            break;
            
        case 8:
            signualar_quad_weight_file_ = var_val; 
            break;
            
        case 9:
            regular_quad_weigh_file_ = var_val;
            break;   
            
        case 10:
            w_sph_file_ = var_val; break;                
    
        default:
            cout << "'" << var_name
                 << "' is an invalid parameter."
                 << endl;
            
            break;
    }
}

ostream& operator<<(ostream& output, const FileList& par)
{
    
    output<<"\n ----------------------------------------------"<<endl;
    output<<"  File list"<<endl;
    output<<" ----------------------------------------------"<<endl;
    output<<"  initial_shape_file         : "<<par.initial_shape_file_<<endl;                 
    output<<"  simulation_out_file        : "<<par.simulation_out_file_<<endl;        
    output<<"  leg_transform_file         : "<<par.leg_transform_file_<<endl;                 
    output<<"  inv_leg_transform_file     : "<<par.inv_leg_transform_file_<<endl;         
    output<<"  inv_d1_leg_transform_file  : "<<par.inv_d1_leg_transform_file_<<endl;   
    output<<"  inv_d2_leg_transform_file  : "<<par.inv_d2_leg_transform_file_<<endl;   
    output<<"  all_rotation_matrices_file : "<<par.all_rotation_matrices_file_<<endl; 
    output<<"  signualar_quad_weight_file : "<<par.signualar_quad_weight_file_<<endl; 
    output<<"  regular_quad_weigh_file    : "<<par.regular_quad_weigh_file_<<endl;       
    output<<"  w_sph_file                 : "<<par.w_sph_file_<<endl;                                 
    output<<" ----------------------------------------------"<<endl<<endl;
    
    return output;
}    

template <typename T>
struct AllParams 
{
    SurfaceParams<T> s_par;
    TimeStepperParams<T> t_par;
    FileList f_list;
    
    AllParams();
    void SetMember(string var_name, string var_val);

  private:
    map<string, int> mapStringValues;
};

template<typename T>
AllParams<T>::AllParams()
{
    this->mapStringValues["SurfaceParams"]     = 1;
    this->mapStringValues["TimeStepperParams"] = 2;
    this->mapStringValues["FileList"]          = 3;
}

template<typename T>
void AllParams<T>::SetMember(string var_name, string var_val)
{
    string::size_type i = var_name.find_first_of ( "." );
    string line;
    if ( i != string::npos)
    {
        line = var_name.substr(0,i);
        var_name = var_name.substr(++i);
    }
    else
    {
        line = var_name;
    }
    
    switch (this->mapStringValues[line])
    {
        case 1:
            this->s_par.SetMember(var_name,var_val);
            break;
            
        case 2:
            this->t_par.SetMember(var_name,var_val);
            break;

        case 3:
            this->f_list.SetMember(var_name,var_val);
            break;

        default:
            cout << "'" << var_name
                 << "' is an invalid parameter."
                 << endl;
            
            break;
    }
}

#endif //_TIMESTEPPER_H_



