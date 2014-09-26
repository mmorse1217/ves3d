template <typename Container>
OperatorsMats<Container>::OperatorsMats(bool readFromFile,
    const Parameters<value_type> &params) :
    p_(params.sh_order),
    p_up_(params.rep_up_freq),
    data_(getDataLength(params)),
    mats_p_(Container::getDevice(), p_, data_.begin(), readFromFile),
    mats_p_up_(Container::getDevice(), p_up_, data_.begin() +
        SHTMats<value_type, device_type>::getDataLength(p_), readFromFile)
{
    int np = 2 * p_ * ( p_ + 1);
    int np_up = 2 * p_up_ * (p_up_ + 1);

    int rot_mat_size =  0;
    int spharm_rot_size = 0;

    if (params.singular_stokes == Direct ||
        params.singular_stokes == DirectEagerEval)
        rot_mat_size =  np * np * (p_ + 1);
    else
        spharm_rot_size = (p_ + 1) * (p_ * (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);

    quad_weights_       = data_.begin() +
        SHTMats<value_type,device_type>::getDataLength(p_) +
        SHTMats<value_type,device_type>::getDataLength(p_up_);

    quad_weights_p_up_  = quad_weights_       + np;
    sing_quad_weights_  = quad_weights_p_up_  + np_up;
    w_sph_              = sing_quad_weights_  + np;
    all_rot_mats_       = w_sph_              + np;
    sh_rot_mats_        = all_rot_mats_       + rot_mat_size;

    if (rot_mat_size == 0)
        all_rot_mats_ = NULL;
    else
        sh_rot_mats_ = NULL;

    DataIO fileIO;

    if(readFromFile)
    {
        std::string tname;
        char fname[500];

        if ( typeid(value_type) == typeid(float) )
            tname = "single";
        else
            tname = "double";

        sprintf(fname,"precomputed/quad_weights_%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, quad_weights_ - data_.begin(), np);

        sprintf(fname,"precomputed/sing_quad_weights_%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, sing_quad_weights_- data_.begin() , np);

        sprintf(fname,"precomputed/w_sph_%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, w_sph_- data_.begin(), np);

        if (rot_mat_size != 0)
        {
            sprintf(fname,"precomputed/all_rot_mats_%u_%s.txt",p_,tname.c_str());
            fileIO.ReadData(fname, data_, DataIO::ASCII, all_rot_mats_- data_.begin(),
                rot_mat_size);
        } else {
            sprintf(fname,"precomputed/SpHarmRotMats_p%u_%s.txt",p_,tname.c_str());
            fileIO.ReadData(fname, data_, DataIO::ASCII, sh_rot_mats_- data_.begin(),
                spharm_rot_size);
        }

        //p
        sprintf(fname,"precomputed/legTrans%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_.dlt_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(fname,"precomputed/legTransInv%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_.dlt_inv_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(fname,"precomputed/d1legTrans%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_.dlt_inv_d1_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(fname,"precomputed/d2legTrans%u_%s.txt",p_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_.dlt_inv_d2_- data_.begin(),
            mats_p_.getDLTLength());

        //p_up
        sprintf(fname,"precomputed/quad_weights_%u_%s.txt",p_up_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, quad_weights_p_up_ - data_.begin(), np_up);

        sprintf(fname,"precomputed/legTrans%u_%s.txt",p_up_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_up_.dlt_- data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(fname,"precomputed/legTransInv%u_%s.txt",p_up_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_up_.dlt_inv_ - data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(fname,"precomputed/d1legTrans%u_%s.txt",p_up_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_up_.dlt_inv_d1_ - data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(fname,"precomputed/d2legTrans%u_%s.txt",p_up_,tname.c_str());
        fileIO.ReadData(fname, data_, DataIO::ASCII, mats_p_up_.dlt_inv_d2_ - data_.begin(),
            mats_p_up_.getDLTLength());
    }
}

template <typename Container>
size_t OperatorsMats<Container>::getDataLength(const Parameters<value_type>
    &params) const
{
    int np = 2 * p_ * ( p_ + 1);
    int np_up = 2 * p_up_ * ( p_up_ + 1);
    int rot_mat_size = 0;
    int spharm_rot_size = 0;

    if (params.singular_stokes == Direct ||
        params.singular_stokes == DirectEagerEval)
        rot_mat_size =  np * np * (p_ + 1);
    else
        spharm_rot_size = (p_ + 1) * (p_ *
            (4 * p_ * p_ -  1)/3 + 4 * p_ * p_);

    return(3*np + np_up + rot_mat_size +  spharm_rot_size +
        SHTMats<value_type, device_type>::getDataLength(p_)      +
        SHTMats<value_type, device_type>::getDataLength(p_up_));
}
