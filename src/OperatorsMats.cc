template <typename T, typename IO>
OperatorsMats<T, IO>::OperatorsMats(IO &fileIO_in, bool readFromFile, 
    const Parameters<T> &params) :
    fileIO_(fileIO_in),
    p_(params.sh_order), 
    p_up_(params.rep_up_freq),
    data_((T*) fileIO_.device_.Malloc(getDataLength() * sizeof(T))),
    mats_p_(p_, data_, readFromFile),
    mats_p_up_(p_up_, data_ + SHTMats<T>::getDataLength(p_), readFromFile)
{
    int np = 2 * p_ * ( p_ + 1);   
    
    quad_weights_       = data_ + 
        SHTMats<T>::getDataLength(p_) + 
        SHTMats<T>::getDataLength(p_up_);
    
    sing_quad_weights_  = quad_weights_       + np;
    w_sph_              = sing_quad_weights_  + np;
    all_rot_mats_       = w_sph_              + np;   

    int rot_mat_size =  np * np * (p_ + 1);
    
    if(readFromFile)
    {
        char fname[500];

        sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, quad_weights_);

        sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, sing_quad_weights_);

        sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        fileIO_.ReadData(fname, np, w_sph_);

        sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p_);
        fileIO_.ReadData(fname, rot_mat_size, all_rot_mats_);

        //p
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_.getDLTLength(), mats_p_.dlt_);

        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_.getDLTLength(), mats_p_.dlt_inv_);

        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_.getDLTLength(), mats_p_.dlt_inv_d1_);
        
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_.getDLTLength(), mats_p_.dlt_inv_d2_);
        
        //p_up
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_);

        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_);

        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_d1_);
        
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_d2_);
    }
}

template <typename T, typename IO>
OperatorsMats<T, IO>::~OperatorsMats()
{
    fileIO_.device_.Free(data_);
}

template <typename T, typename IO>
size_t OperatorsMats<T, IO>::getDataLength()
{
    int np = 2 * p_ * ( p_ + 1);
    int rot_mat_size =  np * np * (p_ + 1);
    
    return(3*np + rot_mat_size +  SHTMats<T>::getDataLength(p_) + 
        SHTMats<T>::getDataLength(p_up_));
}
