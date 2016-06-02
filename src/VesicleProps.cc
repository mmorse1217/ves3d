template<typename T>
VesicleProperties<T>::VesicleProperties() :
    has_contrast(false),
    has_excess_density(false)
{
    set_name("vesicle_properties");
    bending_modulus    .set_name("bending_modulus");
    viscosity_contrast .set_name("viscosity_contrast");
    excess_density     .set_name("excess_density");
    vel_coeff          .set_name("vel_coeff");
    dl_coeff           .set_name("dl_coeff");
    bending_coeff      .set_name("bending_coeff");
}

template<typename T>
int VesicleProperties<T>::n_props(3);

template<typename T>
Error_t VesicleProperties<T>::update(){

    size_type nves(viscosity_contrast.size());
    vel_coeff     .resize(nves);
    dl_coeff      .resize(nves);
    bending_coeff .resize(nves);

    value_type *buffer = new value_type[nves];

    // setting the double layer and velocity coefficients
    for (size_type iV(0); iV<nves; ++iV) buffer[iV] = 0.5;
    vel_coeff.getDevice().Memcpy(vel_coeff.begin(),buffer,
        nves*sizeof(value_type),device_type::MemcpyHostToDevice);
    vel_coeff.getDevice().axpy(0.5, viscosity_contrast.begin(),
        vel_coeff.begin(), nves, vel_coeff.begin());

    for (size_type iV(0); iV<nves; ++iV) buffer[iV] = 1.0;
    dl_coeff.getDevice().Memcpy(dl_coeff.begin(),buffer,
        nves*sizeof(value_type),device_type::MemcpyHostToDevice);
    dl_coeff.getDevice().axpy(-1.0, viscosity_contrast.begin(),
        dl_coeff.begin(), nves, dl_coeff.begin());

    // bending coefficient
    bending_modulus.getDevice().axpy<value_type>(-1.0,bending_modulus.begin(),
        NULL, nves, bending_coeff.begin());

    // check contrast and excess density to set flags
    has_contrast=false;
    viscosity_contrast.getDevice().Memcpy(buffer,viscosity_contrast.begin(),
        nves*sizeof(value_type),device_type::MemcpyDeviceToHost);
    for (size_type iV(0); iV<nves; ++iV) has_contrast |= (fabs(buffer[iV]-1.0)>1e-12);

    has_excess_density=false;
    excess_density.getDevice().Memcpy(buffer,excess_density.begin(),
        nves*sizeof(value_type),device_type::MemcpyDeviceToHost);
    for (size_type iV(0); iV<nves; ++iV) has_excess_density |= (fabs(buffer[iV])>1e-12);

    delete[] buffer;

    return ErrorEvent::Success;
}

template<typename T>
T* VesicleProperties<T>::getPropIdx(int i){
    ASSERT(i<n_props,"Index out of bound");

    switch (i) {
        case 0:
            return &bending_modulus;
        case 1:
            return &viscosity_contrast;
        default: /* case 2: */
            return &excess_density;
    }
}

template<typename T>
Error_t VesicleProperties<T>::setFromParams(Parameters<value_type> &params){

    INFO("Populating vesicle properties from parameters");
    int nves(params.n_surfs);
    value_type *buffer = new value_type[nves];

    bending_modulus.resize(nves);
    for (int i(0);i<nves;++i) buffer[i]=params.bending_modulus;
    bending_modulus.getDevice().Memcpy(bending_modulus.begin(),buffer,
        nves * sizeof(value_type), T::device_type::MemcpyHostToDevice);

    excess_density.resize(nves);
    for (int i(0);i<nves;++i) buffer[i]=params.excess_density;
    excess_density.getDevice().Memcpy(excess_density.begin(),buffer,
        nves * sizeof(value_type),T::device_type::MemcpyHostToDevice);

    viscosity_contrast.resize(nves);
    for (int i(0);i<nves;++i) buffer[i]=params.viscosity_contrast;
    viscosity_contrast.getDevice().Memcpy(viscosity_contrast.begin(),buffer,
        nves * sizeof(value_type), T::device_type::MemcpyDeviceToDevice);

    delete buffer;
    // check
    for (int iP(0);iP<n_props;++iP){
        T* prp(getPropIdx(iP));
        ASSERT(prp->size()==nves,"property has wrong size (maybe uninitialized)");
    }

    update();

    return ErrorEvent::Success;
}

template<typename T>
Error_t VesicleProperties<T>::pack(std::ostream &os, Format format) const
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    os<<"VESICLEPROPS\n";
    os<<"version: "<<VERSION<<"\n";
    os<<"name: "<<Streamable::name_<<"\n";
    bending_modulus   .pack(os,format);
    viscosity_contrast.pack(os,format);
    excess_density    .pack(os,format);
    os<<"/VESICLEPROPS\n";
    return ErrorEvent::Success;
}

template<typename T>
Error_t VesicleProperties<T>::unpack(std::istream &is, Format format)
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s,key;
    int version(0);

    is>>s;
    ASSERT(s=="VESICLEPROPS", "Bad input string (missing header).");

    is>>key;
    if (key=="version:") {
        is>>version>>key;
        if (key=="+") {++version;is>>key;}
    };
    is>>Streamable::name_;
    ASSERT(key=="name:", "bad key name");

    bending_modulus   .unpack(is,format);
    viscosity_contrast.unpack(is,format);
    excess_density    .unpack(is,format);

    is>>s;
    ASSERT(s=="/VESICLEPROPS", "Bad input string (missing footer).");
    update();
    INFO("Unpacked "<<Streamable::name_<<" data from version "<<version<<" (current version "<<VERSION<<")");

    return ErrorEvent::Success;
}
