template <typename Container>
OperatorsMats<Container>::OperatorsMats(bool readFromFile, 
    const Parameters<value_type> &params) :
    p_(params.sh_order), 
    p_up_(params.rep_up_freq),
    data_(getDataLength()),
    mats_p_(Container::getDevice(), p_, data_.begin(), readFromFile),
    mats_p_up_(Container::getDevice(), p_up_, data_.begin() + 
        SHTMats<value_type, device_type>::getDataLength(p_), readFromFile)
{
    int np = 2 * p_ * ( p_ + 1);   
    int rot_mat_size =  np * np * (p_ + 1);
    int spharm_rot_size = (p_ + 1) * (p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);
    
    quad_weights_       = data_.begin() + 
        SHTMats<value_type,device_type>::getDataLength(p_) + 
        SHTMats<value_type,device_type>::getDataLength(p_up_);
    
    sing_quad_weights_  = quad_weights_       + np;
    w_sph_              = sing_quad_weights_  + np;
    all_rot_mats_       = w_sph_              + np;   
    sh_rot_mats_        = all_rot_mats_       + rot_mat_size;

    DataIO fileIO;

    if(readFromFile)
    {
        char fname[500];
        sprintf(fname,"precomputed/quad_weights_%u_single.txt",p_);
        fileIO.ReadData(fname, data_,quad_weights_ - data_.begin(), np);

        sprintf(fname,"precomputed/sing_quad_weights_%u_single.txt",p_);
        fileIO.ReadData(fname, data_, sing_quad_weights_- data_.begin() , np);

        sprintf(fname,"precomputed/w_sph_%u_single.txt",p_);
        fileIO.ReadData(fname, data_, w_sph_- data_.begin(), np);

        sprintf(fname,"precomputed/all_rot_mats_%u_single.txt",p_);
        fileIO.ReadData(fname, data_, all_rot_mats_- data_.begin(),
            rot_mat_size);
        
//         sprintf(fname,"precomputed/SpHarmRotMats_p%u_float.txt",p_);
//         fileIO.ReadData(fname, data_, sh_rot_mats_- data_.begin(),
//             spharm_rot_size);

        //p
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_);
        fileIO.ReadData(fname, data_, mats_p_.dlt_ - data_.begin(), 
            mats_p_.getDLTLength());

        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_);
        fileIO.ReadData(fname, data_, mats_p_.dlt_inv_ - data_.begin(), 
            mats_p_.getDLTLength());

        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_);
        fileIO.ReadData(fname, data_, mats_p_.dlt_inv_d1_ - data_.begin(),
            mats_p_.getDLTLength());
        
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_);
        fileIO.ReadData(fname, data_, mats_p_.dlt_inv_d2_- data_.begin(),
            mats_p_.getDLTLength());
        
        //p_up
        sprintf(fname,"precomputed/legTrans%u_single.txt",p_up_);
        fileIO.ReadData(fname, data_, mats_p_up_.dlt_- data_.begin(), 
            mats_p_up_.getDLTLength());
        
        sprintf(fname,"precomputed/legTransInv%u_single.txt",p_up_);
        fileIO.ReadData(fname, data_, mats_p_up_.dlt_inv_ - data_.begin(), 
            mats_p_up_.getDLTLength());

        sprintf(fname,"precomputed/d1legTrans%u_single.txt",p_up_);
        fileIO.ReadData(fname, data_, mats_p_up_.dlt_inv_d1_ - data_.begin(),
            mats_p_up_.getDLTLength());
        
        sprintf(fname,"precomputed/d2legTrans%u_single.txt",p_up_);
        fileIO.ReadData(fname, data_, mats_p_up_.dlt_inv_d2_ - data_.begin(),
            mats_p_up_.getDLTLength());
    }
}

template <typename Container>
size_t OperatorsMats<Container>::getDataLength() const
{
    int np = 2 * p_ * ( p_ + 1);
    int rot_mat_size =  np * np * (p_ + 1);
    int spharm_rot_size = (p_ + 1) * (p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);
    
    return(3*np + rot_mat_size +  spharm_rot_size + 
        SHTMats<value_type, device_type>::getDataLength(p_)      + 
        SHTMats<value_type, device_type>::getDataLength(p_up_));
}
