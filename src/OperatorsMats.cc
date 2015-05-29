template <typename Container>
OperatorsMats<Container>::OperatorsMats(bool readFromFile,
    const Parameters<value_type> &params) :
    p_(params.sh_order),
    p_up_(params.upsample_freq),
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

    value_type* data_ptr = data_.begin() +
        SHTMats<value_type,device_type>::getDataLength(p_) +
        SHTMats<value_type,device_type>::getDataLength(p_up_);

    quad_weights_         =data_ptr; data_ptr+= np;
    quad_weights_p_up_    =data_ptr; data_ptr+= np_up;
    sing_quad_weights_    =data_ptr; data_ptr+= np;
    sing_quad_weights_up_ =data_ptr; data_ptr+= np_up;
    w_sph_                =data_ptr; data_ptr+= np;
    w_sph_up_             =data_ptr; data_ptr+= np_up;
    all_rot_mats_         =data_ptr; data_ptr+= rot_mat_size;
    sh_rot_mats_          =data_ptr; data_ptr+= spharm_rot_size;
    assert((data_ptr-data_.begin())==getDataLength(params));

    if (rot_mat_size == 0)
        all_rot_mats_ = NULL;
    else
        sh_rot_mats_ = NULL;

    DataIO fileIO;

    if(readFromFile)
    {
        std::string tname;
        char buffer[500];

        if ( typeid(value_type) == typeid(float) )
            tname = "single";
        else
            tname = "double";

        sprintf(buffer,"precomputed/quad_weights_%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, quad_weights_ - data_.begin(), np);

        sprintf(buffer,"precomputed/sing_quad_weights_%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, sing_quad_weights_- data_.begin() , np);

        sprintf(buffer,"precomputed/w_sph_%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, w_sph_- data_.begin(), np);

        if (rot_mat_size != 0)
        {
            sprintf(buffer,"precomputed/all_rot_mats_%u_%s.txt", p_,tname.c_str());
            fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, all_rot_mats_- data_.begin(),
                rot_mat_size);
        } else {
            sprintf(buffer,"precomputed/SpHarmRotMats_p%u_%s.txt", p_,tname.c_str());
            fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, sh_rot_mats_- data_.begin(),
                spharm_rot_size);
        }

        //p
        sprintf(buffer,"precomputed/legTrans%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_.dlt_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(buffer,"precomputed/legTransInv%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_.dlt_inv_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(buffer,"precomputed/d1legTrans%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_.dlt_inv_d1_ - data_.begin(),
            mats_p_.getDLTLength());

        sprintf(buffer,"precomputed/d2legTrans%u_%s.txt", p_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_.dlt_inv_d2_- data_.begin(),
            mats_p_.getDLTLength());

        //p_up
        sprintf(buffer,"precomputed/quad_weights_%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, quad_weights_p_up_ - data_.begin(), np_up);

        sprintf(buffer,"precomputed/sing_quad_weights_%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, sing_quad_weights_up_- data_.begin() , np_up);

        sprintf(buffer,"precomputed/w_sph_%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, w_sph_up_- data_.begin(), np_up);

        sprintf(buffer,"precomputed/legTrans%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_up_.dlt_- data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(buffer,"precomputed/legTransInv%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_up_.dlt_inv_ - data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(buffer,"precomputed/d1legTrans%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_up_.dlt_inv_d1_ - data_.begin(),
            mats_p_up_.getDLTLength());

        sprintf(buffer,"precomputed/d2legTrans%u_%s.txt", p_up_,tname.c_str());
        fileIO.ReadData(FullPath(buffer), data_, DataIO::ASCII, mats_p_up_.dlt_inv_d2_ - data_.begin(),
            mats_p_up_.getDLTLength());

	INFO("Matrices are loaded from file");
    } else {
      INFO("Object created with no data");
    }
}

template <typename Container>
size_t OperatorsMats<Container>::getDataLength(const Parameters<value_type> &params) const
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

    return(3*np + 3*np_up + rot_mat_size +  spharm_rot_size +
        SHTMats<value_type, device_type>::getDataLength(p_)      +
        SHTMats<value_type, device_type>::getDataLength(p_up_));
}

template <typename Container>
const SHTMats<typename Container::value_type, typename Container::device_type>&
OperatorsMats<Container>::getShMats(int order) const
{
    if (order==p_)
	return mats_p_;
    else /* if (order==p_up_) */
	return mats_p_up_;
}
