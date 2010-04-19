#ifndef _TIMESTEPPER_H_
#define _TIMESTEPPER_H_

#include "Surface.h"
#include "VesUtil.h"
#include "DataIO.h"
//#include "Monitor.h"

template<typename T> 
class TimeStepperMult
{
  public:
    // {{{ Member variables
    T ts_;
    bool saveData;
    bool verbose;
    int n_steps_;
    int num_ves_;
    DataIO<T> &fileIO;
    Surface<T> **vesicle_vec;
    Vectors<T> **velocity_, **vel_bending_, **vel_tension_;
    void *user;

    VelField<T> &bg_flow_;
    //Monitor<T> &mntr;
    void(*Interaction_)(Surface<T> **, int, void*);
    bool(*userMonitor)(TimeStepperMult<T>&, int current_time_step);

    T **quad_weights_;
    // }}}
    
    // {{{ Constructor 
    TimeStepperMult(T ts_in, int n_steps_in, Surface<T> **ves_in, int num_ves_in, 
        DataIO<T> &fileIO_in, VelField<T> &bg_flow_in, 
        T* quad_weights_p_up_in,
        void(*Interaction_in)(Surface<T> **, int, void*)) : 
        ts_(ts_in),
        n_steps_(n_steps_in), 
        fileIO(fileIO_in),
        saveData(true),
        verbose(true),
        userMonitor(NULL),
        vesicle_vec(ves_in),
        num_ves_(num_ves_in),
        bg_flow_(bg_flow_in),
        //mntr(mntr_in),
        velocity_(new Vectors<T>*[num_ves_]),
        vel_bending_(new Vectors<T>*[num_ves_]),
        vel_tension_(new Vectors<T>*[num_ves_]),
        Interaction_(Interaction_in),
        quad_weights_(new T*[num_ves_])
    {
        for(int ii=0;ii<num_ves_;++ii)
        {
            velocity_[ii]    = new Vectors<T>(vesicle_vec[ii]->device_, vesicle_vec[ii]->params_.p_, vesicle_vec[ii]->params_.n_surfs_);
            vel_bending_[ii] = new Vectors<T>(vesicle_vec[ii]->device_, vesicle_vec[ii]->params_.p_, vesicle_vec[ii]->params_.n_surfs_);
            vel_tension_[ii] = new Vectors<T>(vesicle_vec[ii]->device_, vesicle_vec[ii]->params_.p_, vesicle_vec[ii]->params_.n_surfs_);
            
            int np_up = 2 * vesicle_vec[ii]->params_.rep_up_freq_ * (vesicle_vec[ii]->params_.rep_up_freq_ + 1);
            quad_weights_[ii] = vesicle_vec[ii]->device_.Malloc(np_up);
            vesicle_vec[ii]->device_.Memcpy(quad_weights_[ii], quad_weights_p_up_in, np_up, MemcpyHostToDevice);
        }
    };
    // }}}

    // {{{ Destructor
    ~TimeStepperMult()
    {
        for(int ii=0;ii<num_ves_;++ii)
        {
            delete velocity_[ii];
            delete vel_bending_[ii];
            delete vel_tension_[ii];
            vesicle_vec[ii]->device_.Free(quad_weights_[ii]);
        }
        
        delete[] velocity_;
        delete[] vel_bending_;
        delete[] vel_tension_;
        delete[] quad_weights_;
    }
    // }}}

    // {{{ Resize
    ///@todo Resize is obsolete 
    void Resize(int new_size)
    {
        this->vesicle_.Resize(new_size);
        this->velocity_.Resize(new_size);
        this->vel_bending_.Resize(new_size);
        this->vel_tension_.Resize(new_size);
    }
    
    // }}}

    // {{{ Evolve in time
    void EvolveInTime(int n){ n_steps_=n; 	EvolveInTime();} 

    void EvolveInTime()
    {
        if(verbose){
            for(int ii=0;ii<num_ves_;++ii)
                cout<<vesicle_vec[ii]->params_;
            
            cout<<" ------------------------------------"<<endl;
            cout<<"  Time stepper parameters"<<endl;
            cout<<" ------------------------------------"<<endl;
            cout<<"  ts                : "<<ts_<<endl;
            cout<<"  n_steps           : "<<n_steps_<<endl;
            cout<<" ------------------------------------"<<endl<<endl;
        }

        for(int current_time_step=0;current_time_step<n_steps_;current_time_step++)
        {

#pragma omp parallel num_threads(num_ves_)
            {
                assert(omp_get_num_threads() == num_ves_);
                int this_thread = omp_get_thread_num();
                Surface<T> *this_vesicle(vesicle_vec[this_thread]);
                
#pragma omp critical
                if(verbose){ cout<<"t = "<<current_time_step<<", vesicle #"<<this_thread<<endl; }
                
                //Set up the current state
                this_vesicle->UpdateAll();
                
                // {{{ Interaction
                if(Interaction_ !=NULL)
                {
                    //Interaction
                    avpw(this_vesicle->tension_, this_vesicle->tensile_force_, this_vesicle->bending_force_, *(vel_tension_[this_thread]));
                    
                    //vel_tension_ holds the interaction force
                    xvpb(this_vesicle->w_,*(vel_tension_[this_thread]), (T) 0.0, *(vel_tension_[this_thread]));
                
                    //up-sampling x_ to V10
                    this_vesicle->device_.InterpSh(this_vesicle->params_.p_, 3*this_vesicle->params_.n_surfs_,
                         this_vesicle->x_.data_, this_vesicle->work_arr, this_vesicle->shc, this_vesicle->params_.rep_up_freq_, this_vesicle->V10.data_);
                
                    //up-sampling vel_tension to V11
                    this_vesicle->device_.InterpSh(this_vesicle->params_.p_, 3*this_vesicle->params_.n_surfs_,
                        vel_tension_[this_thread]->data_, this_vesicle->work_arr, this_vesicle->shc, this_vesicle->params_.rep_up_freq_, this_vesicle->V11.data_);

                    //the self-interaction -- V12
                    this_vesicle->device_.DirectStokes(this_vesicle->V10.GetFunLength(), this_vesicle->params_.n_surfs_, 
                        0, this_vesicle->V10.GetFunLength(), quad_weights_[this_thread], this_vesicle->V10.data_,
                         this_vesicle->V10.data_, this_vesicle->V11.data_, this_vesicle->V12.data_);
            
                    ///@todo the multiplication by the quadrature weights is slow
                    int stride = this_vesicle->V11.GetFunLength();
                    for(int ii=0;ii<this_vesicle->params_.n_surfs_;++ii)
                    {
                         this_vesicle->device_.xvpb(quad_weights_[this_thread], this_vesicle->V11.data_ + 3*ii*stride,
                             (T) 0.0, stride, 1, this_vesicle->V11.data_ + 3*ii*stride);
                     }
                
                     // Call to fast summation algorithm for this_vesicle--this_vesicle interactions. 
                     //Interaction to V13 
#pragma omp barrier
#pragma omp master
                    {
                        Interaction_(vesicle_vec, num_ves_, user);
                    }
#pragma omp barrier

                     axpy((T) -1.0, this_vesicle->V12, this_vesicle->V13, this_vesicle->V13);

                     //filtering to 1/3
                     this_vesicle->device_.ShAna(this_vesicle->V13.data_, this_vesicle->work_arr, this_vesicle->params_.rep_up_freq_, 
                         3*this_vesicle->params_.n_surfs_, this_vesicle->V11.data_);
            
                     this_vesicle->device_.ScaleFreqs(this_vesicle->params_.rep_up_freq_, 3*this_vesicle->params_.n_surfs_, 
                         this_vesicle->V11.data_, this_vesicle->alpha_q, this_vesicle->V11.data_);

                     this_vesicle->device_.Resample(this_vesicle->params_.rep_up_freq_, 3*this_vesicle->params_.n_surfs_, 
                         this_vesicle->params_.p_, this_vesicle->V11.data_, this_vesicle->shc);

                     this_vesicle->device_.ShSyn(this_vesicle->shc, this_vesicle->work_arr, this_vesicle->params_.p_, 
                         3*this_vesicle->params_.n_surfs_, velocity_[this_thread]->data_);
                
                     //Background flow
                     bg_flow_.GetVel(this_vesicle->x_, this_vesicle->S1, *(vel_tension_[this_thread]));
                     axpy((T) 1.0, *(vel_tension_[this_thread]), *(velocity_[this_thread]), *(velocity_[this_thread]));
                } 
                else
                {     
                    bg_flow_.GetVel(this_vesicle->x_, this_vesicle->S1, *(velocity_[this_thread]));
                }
                // }}}
                    
                //Calculate stokes
                this_vesicle->StokesMatVec(this_vesicle->bending_force_, *(vel_bending_[this_thread]));
                axpy((T) 1.0, this_vesicle->bending_force_, *(velocity_[this_thread]), *(velocity_[this_thread]));
                this_vesicle->StokesMatVec(this_vesicle->tensile_force_, *(vel_tension_[this_thread]));
                    
                //Calculate tension
                this_vesicle->GetTension(*(velocity_[this_thread]), *(vel_tension_[this_thread]), this_vesicle->tension_);
                this_vesicle->device_.avpw(this_vesicle->tension_, vel_tension_[this_thread]->data_,
                    velocity_[this_thread]->data_, velocity_[this_thread]->GetFunLength(), this_vesicle->params_.n_surfs_, velocity_[this_thread]->data_);

                //Advance in time
                axpy(ts_, *(velocity_[this_thread]), this_vesicle->x_, this_vesicle->x_);
            
                this_vesicle->Reparam();
            }
            if ( userMonitor!=NULL )
            { 
                bool flag = userMonitor(*this,current_time_step); 
                if(flag) break;
            }
        }
    }
    // }}}
};


//// Parameters
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
