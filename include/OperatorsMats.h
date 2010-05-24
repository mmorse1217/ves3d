#ifndef _OPERATORSMATS_H_
#define _OPERATORSMATS_H_

#include "DataIO.h"

template <typename T>
struct OperatorsMats
{
    DataIO<T,CPU> &fileIO_;
    int p_;
    int p_up_;
    
    T *data_; // holds all arrays one after another
    T *quad_weights_;
    T *all_rot_mats_;
    T *sing_quad_weights_;
    T *w_sph_;
    
    T *leg_trans_p_;
    T *leg_trans_inv_p_;
    T *d1_leg_trans_p_;
    T *d2_leg_trans_p_;

    T *quad_weights_p_up_;
    T *leg_trans_p_up_;
    T *leg_trans_inv_p_up_;
    T *d1_leg_trans_p_up_;
    T *d2_leg_trans_p_up_;

    long int GetDataLength();
    OperatorsMats(DataIO<T,CPU> &fileIO_in, int p_in, int p_up_in, bool readFromFile);
    ~OperatorsMats();

  private:
    OperatorsMats(const OperatorsMats<T>& mat_in);
    OperatorsMats<T>& operator=(const OperatorsMats<T>& vec_in);
};

template <typename T>
long int OperatorsMats<T>::GetDataLength()
{
    int np = 2 * p_ * ( p_ + 1);
    int rot_mat_size =  np * np * (p_ + 1);
    int np_up = 2 * p_up_ * ( p_up_ + 1);
    int leg_size = (p_ + 1) * (p_+1) * (p_ +2);
    int leg_size_up = (p_up_ + 1) * (p_up_+1) * (p_up_ +2);

    return(3*np + 4*leg_size + 4*leg_size_up + rot_mat_size + np_up);
}

template <typename T>
OperatorsMats<T>::OperatorsMats(DataIO<T,CPU> &fileIO_in, int p_in, int p_up_in, bool readFromFile) :
    fileIO_(fileIO_in),
    p_(p_in), 
    p_up_(p_up_in)
{
    int np = 2 * p_ * ( p_ + 1);
    int rot_mat_size =  np * np * (p_ + 1);
    int np_up = 2 * p_up_ * ( p_up_ + 1);
    int leg_size = (p_ + 1) * (p_ + 1) * (p_ + 2);
    int leg_size_up = (p_up_ + 1) * (p_up_ + 1) * (p_up_ + 2);
    
    data_ = (T*) fileIO_.device_.Malloc(GetDataLength() * sizeof(T));
    
    quad_weights_       = data_;
    all_rot_mats_       = quad_weights_       + np;
    sing_quad_weights_  = all_rot_mats_       + rot_mat_size;
    w_sph_              = sing_quad_weights_  + np;
    
    leg_trans_p_        = w_sph_              + np;
    leg_trans_inv_p_    = leg_trans_p_        + leg_size;
    d1_leg_trans_p_     = leg_trans_inv_p_    + leg_size;
    d2_leg_trans_p_     = d1_leg_trans_p_     + leg_size;
     
    //p_up_
    quad_weights_p_up_  = d2_leg_trans_p_     + leg_size;
    leg_trans_p_up_     = quad_weights_p_up_  + np_up;
    leg_trans_inv_p_up_ = leg_trans_p_up_     + leg_size_up;
    d1_leg_trans_p_up_  = leg_trans_inv_p_up_ + leg_size_up;
    d2_leg_trans_p_up_  = d1_leg_trans_p_up_  + leg_size_up;

    
    if(readFromFile)
    {
        char fname[500];

        //sprintf(fname,"%s/precomputed/quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_);
        sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, quad_weights_);

        //sprintf(fname,"%s/precomputed/all_rot_mats_%u_single.txt",getenv("VES3D_DIR"),p_);
        sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p_);
        fileIO_.ReadData(fname, rot_mat_size, all_rot_mats_);

        //sprintf(fname,"%s/precomputed/sing_quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_);
        sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, sing_quad_weights_);
    
        //sprintf(fname,"%s/precomputed/w_sph_%u_single.txt",getenv("VES3D_DIR"),p_);
        sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, w_sph_);

//         //sprintf(fname,"%s/precomputed/legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
//         sprintf(fname,"precomputed/legTrans%u_single.txt",p_);
//         fileIO_.ReadData(fname, leg_size, leg_trans_p_);
        
//         //sprintf(fname,"%s/precomputed/legTransInv%u_single.txt",getenv("VES3D_DIR"),p_);
//         sprintf(fname,"precomputed/legTransInv%u_single.txt",p_);
//         fileIO_.ReadData(fname, leg_size, leg_trans_inv_p_);

//         //sprintf(fname,"%s/precomputed/d1legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
//         sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_);
//         fileIO_.ReadData(fname, leg_size, d1_leg_trans_p_);

//         //sprintf(fname,"%s/precomputed/d2legTrans%u_single.txt",getenv("VES3D_DIR"),p_);
//         sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_);
//         fileIO_.ReadData(fname, leg_size, d2_leg_trans_p_);

        //p_up_
        //sprintf(fname,"%s/precomputed/quad_weights_%u_single.txt",getenv("VES3D_DIR"),p_up_);
        sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_up_);
        fileIO_.ReadData(fname, np_up, quad_weights_p_up_);
        
        //sprintf(fname,"%s/precomputed/legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, leg_size, leg_trans_p_up_);

        //sprintf(fname,"%s/precomputed/legTransInv%u_single.txt",getenv("VES3D_DIR"),p_up_);
        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_up_);
        fileIO_.ReadData(fname, leg_size, leg_trans_inv_p_up_);

        //sprintf(fname,"%s/precomputed/d1legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, leg_size, d1_leg_trans_p_up_);

        //sprintf(fname,"%s/precomputed/d2legTrans%u_single.txt",getenv("VES3D_DIR"),p_up_);
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, leg_size, d2_leg_trans_p_up_);
    }
}

template <typename T>
OperatorsMats<T>::~OperatorsMats()
{
    fileIO_.device_.Free(data_);
}

#endif //_OPERATORSMATS_H_
