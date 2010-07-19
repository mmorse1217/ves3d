template <typename T, typename Device>
OperatorsMats<T, Device>::OperatorsMats(const Device &dev, DataIO 
    &fileIO_in, bool readFromFile, const Parameters<T> &params) :
    device_(dev),
    fileIO_(fileIO_in),
    p_(params.sh_order), 
    p_up_(params.rep_up_freq),
    data_((T*) device_.Malloc(getDataLength() * sizeof(T))),
    mats_p_(dev, p_, data_, readFromFile),
    mats_p_up_(dev, p_up_, data_ + SHTMats<T, Device>::getDataLength(p_), readFromFile)
{
    int np = 2 * p_ * ( p_ + 1);   
    int rot_mat_size =  np * np * (p_ + 1);
    int spharm_rot_size = (p_ + 1) * (p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);

    quad_weights_       = data_ + 
        SHTMats<T,Device>::getDataLength(p_) + 
        SHTMats<T,Device>::getDataLength(p_up_);
    
    sing_quad_weights_  = quad_weights_       + np;
    w_sph_              = sing_quad_weights_  + np;
    all_rot_mats_       = w_sph_              + np;   
    sh_rot_mats_        = all_rot_mats_ + rot_mat_size;

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
        
        //sprintf(fname,"precomputed/SpHarmRotMats_p%u_float.txt",p_);
        //fileIO_.ReadData(fname, spharm_rot_size, sh_rot_mats_);

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
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_);

        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_up_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_);

        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_d1_);
        
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_up_);
        fileIO_.ReadData(fname, mats_p_up_.getDLTLength(), mats_p_up_.dlt_inv_d2_);
    }
}

template <typename T, typename Device>
OperatorsMats<T, Device>::~OperatorsMats()
{
    device_.Free(data_);
}

template <typename T, typename Device>
size_t OperatorsMats<T, Device>::getDataLength() const
{
    int np = 2 * p_ * ( p_ + 1);
    int rot_mat_size =  np * np * (p_ + 1);
    int spharm_rot_size = (p_ + 1) * (p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);

    return(3*np + rot_mat_size +  spharm_rot_size + 
        SHTMats<T,Device>::getDataLength(p_)      + 
        SHTMats<T,Device>::getDataLength(p_up_));
}

template <typename T, typename Device>
const Device& OperatorsMats<T, Device>::getDevice() const
{
    return(device_);
}
