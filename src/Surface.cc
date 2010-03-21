/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */

template <typename T> 
Surface<T>::Surface(Device<T> &device_in) :
    device_(device_in), 
    x_(device_),
    normal_(device_), 
    h_(device_), 
    w_(device_),
    k_(device_), 
    cu_(device_), 
    cv_(device_),
    bending_force_(device_),
    tensile_force_(device_),
    tension_(0),
    all_rot_mats_(0),
    S1(device_),
    S2(device_),
    S3(device_),
    S4(device_), 
    S5(device_), 
    S6(device_),
    V1(device_),
    V2(device_),
    shc(0),
    work_arr(0),
    alpha_p(0),
    alpha_q(0),
    V10(device_), 
    V11(device_),
    V12(device_), 
    V13(device_),
    S10(device_),
    rot_mat(0),
    w_sph_(device_)
{}

template <typename T> 
Surface<T>::Surface(Device<T> &device_in, SurfaceParams<T> params_in) :
    device_(device_in), 
    params_(params_in),
    x_(device_,params_.p_,params_.n_surfs_),
    normal_(device_, params_.p_,params_.n_surfs_), 
    h_(device_,params_.p_,params_.n_surfs_), 
    w_(device_,params_.p_,params_.n_surfs_),
    k_(device_,params_.p_,params_.n_surfs_), 
    cu_(device_,params_.p_,params_.n_surfs_), 
    cv_(device_,params_.p_,params_.n_surfs_),
    bending_force_(device_,params_.p_,params_.n_surfs_),
    tensile_force_(device_,params_.p_,params_.n_surfs_),
		//GB: the code below is ugly, must be fixed after sc10
    S1(device_,params_.p_,params_.n_surfs_),
    S2(device_,params_.p_,params_.n_surfs_),
    S3(device_,params_.p_,params_.n_surfs_),
    S4(device_,params_.p_,params_.n_surfs_), 
    S5(device_,params_.p_,params_.n_surfs_), 
    S6(device_,params_.p_,params_.n_surfs_),
    V1(device_,params_.p_,params_.n_surfs_),
    V2(device_,params_.p_,params_.n_surfs_),
    V10(device_,params_.rep_up_freq_,params_.n_surfs_), 
    V11(device_,params_.rep_up_freq_,params_.n_surfs_),
    V12(device_,params_.rep_up_freq_,params_.n_surfs_),
    V13(device_,params_.rep_up_freq_,params_.n_surfs_),
    S10(device_,params_.rep_up_freq_,params_.n_surfs_),
    w_sph_(device_,params_.p_,params_.n_surfs_)
{
    shc      = device_.Malloc(6  * params_.rep_up_freq_ *(params_.rep_up_freq_ + 1) * params_.n_surfs_);
    work_arr = device_.Malloc(12 * params_.rep_up_freq_ *(params_.rep_up_freq_ + 1) * params_.n_surfs_);
    alpha_p  = device_.Malloc(params_.p_ *(params_.p_ + 2));
    alpha_q  = device_.Malloc(params_.rep_up_freq_ *(params_.rep_up_freq_ + 2));

    tension_ = device_.Malloc(params_.n_surfs_);
    
    int np = 2 * params_.p_ * (params_.p_ + 1);
    quad_weights_ = device_.Malloc(np);
    all_rot_mats_ = device_.Malloc( np * np * (params_.p_ + 1));
    rot_mat = device_.Malloc(np * np);
    sing_quad_weights_ = device_.Malloc(np);
    vel = device_.Malloc(3*params_.n_surfs_);
    
    T *buffer = (T*) malloc(np * np * (params_.p_ + 1) * sizeof(T));
    //reading quadrature weights and rotation matrix form file
    DataIO<T> myIO;
    char fname[300];
    sprintf(fname,"../data/quad_weights_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(quad_weights_,buffer, np, MemcpyHostToDevice);

    sprintf(fname,"../data/sing_quad_weights_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(sing_quad_weights_, buffer, np, MemcpyHostToDevice);
    
    sprintf(fname,"../data/all_rot_mats_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np * np *(params_.p_ + 1), buffer);
    device_.Memcpy(all_rot_mats_, buffer, np * np *(params_.p_ + 1), MemcpyHostToDevice);
    
    sprintf(fname,"../data/w_sph_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(w_sph_.data_, buffer, np, MemcpyHostToDevice);
    for(int ii=1;ii<params_.n_surfs_;++ii)
        device_.Memcpy(w_sph_.data_ + ii*np, w_sph_.data_, np, MemcpyDeviceToDevice);

    int idx = 0, len;
    for(int ii=0; ii< 2 * params_.p_; ++ii)
    {
        len = params_.p_ + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
            buffer[idx++] = (len-jj)<=(params_.p_-params_.filter_freq_) ? 0 : 1;
    }
    device_.Memcpy(alpha_p, buffer, params_.p_ *(params_.p_ + 2), MemcpyHostToDevice);

    idx = 0;
    for(int ii=0; ii< 2 * params_.rep_up_freq_; ++ii)
    {
        len = params_.rep_up_freq_ + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
            buffer[idx++] = (len-jj)<=(params_.rep_up_freq_-params_.rep_filter_freq_) ? 0 : 1;
    }
    device_.Memcpy(alpha_q, buffer, params_.rep_up_freq_ *(params_.rep_up_freq_ + 2), MemcpyHostToDevice);

    free(buffer);
}

template <typename T> 
Surface<T>::Surface(Device<T> &device_in, SurfaceParams<T> params_in, const Vectors<T> &x_in) :
    device_(device_in), 
    params_(params_in),
    x_(device_,params_.p_,params_.n_surfs_),
    normal_(device_, params_.p_,params_.n_surfs_), 
    h_(device_,params_.p_,params_.n_surfs_), 
    w_(device_,params_.p_,params_.n_surfs_),
    k_(device_,params_.p_,params_.n_surfs_), 
    cu_(device_,params_.p_,params_.n_surfs_), 
    cv_(device_,params_.p_,params_.n_surfs_),
    bending_force_(device_,params_.p_,params_.n_surfs_),
    tensile_force_(device_,params_.p_,params_.n_surfs_),
    S1(device_,params_.p_,params_.n_surfs_),
    S2(device_,params_.p_,params_.n_surfs_),
    S3(device_,params_.p_,params_.n_surfs_),
    S4(device_,params_.p_,params_.n_surfs_), 
    S5(device_,params_.p_,params_.n_surfs_), 
    S6(device_,params_.p_,params_.n_surfs_),
    V1(device_,params_.p_,params_.n_surfs_),
    V2(device_,params_.p_,params_.n_surfs_),
    V10(device_,params_.rep_up_freq_,params_.n_surfs_), 
    V11(device_,params_.rep_up_freq_,params_.n_surfs_),
    V12(device_,params_.rep_up_freq_,params_.n_surfs_),
    V13(device_,params_.rep_up_freq_,params_.n_surfs_),
    S10(device_,params_.rep_up_freq_,params_.n_surfs_),
    w_sph_(device_,params_.p_,params_.n_surfs_)
{
    shc      = device_.Malloc(6  * params_.rep_up_freq_ *(params_.rep_up_freq_ + 1) * params_.n_surfs_);
    work_arr = device_.Malloc(12 * params_.rep_up_freq_ *(params_.rep_up_freq_ + 1) * params_.n_surfs_);
    alpha_p  = device_.Malloc(params_.p_ *(params_.p_ + 2));
    alpha_q  = device_.Malloc(params_.rep_up_freq_ *(params_.rep_up_freq_ + 2));

    tension_ = device_.Malloc(params_.n_surfs_);
    
    int np = 2 * params_.p_ * (params_.p_ + 1);
    quad_weights_ = device_.Malloc(np);
    all_rot_mats_ = device_.Malloc( np * np * (params_.p_ + 1));
    rot_mat = device_.Malloc(np * np);
    sing_quad_weights_ = device_.Malloc(np);
    vel = device_.Malloc(3*params_.n_surfs_);
    
    T *buffer = (T*) malloc(np * np * (params_.p_ + 1) * sizeof(T));
    //reading quadrature weights and rotation matrix form file
    DataIO<T> myIO;
    char fname[300];
    sprintf(fname,"../data/quad_weights_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(buffer, quad_weights_, np, MemcpyHostToDevice);

    sprintf(fname,"../data/sing_quad_weights_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(buffer, sing_quad_weights_, np, MemcpyHostToDevice);
    
    sprintf(fname,"../data/all_rot_mats_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np * np *(params_.p_ + 1), buffer);
    device_.Memcpy(buffer, all_rot_mats_, np * np *(params_.p_ + 1), MemcpyHostToDevice);
    
    sprintf(fname,"../data/w_sph_%u_single.txt",params_.p_);
    myIO.ReadData(fname, np, buffer);
    device_.Memcpy(buffer, w_sph_.data_, np, MemcpyHostToDevice);
    for(int ii=1;ii<params_.n_surfs_;++ii)
        device_.Memcpy(w_sph_.data_ + ii*np, w_sph_.data_, np, MemcpyDeviceToDevice);

    int idx = 0, len;
    for(int ii=0; ii< 2 * params_.p_; ++ii)
    {
        len = params_.p_ + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
            buffer[idx++] = (len-jj)<=(params_.p_-params_.filter_freq_) ? 0 : 1;
    }
    device_.Memcpy(buffer, alpha_p, params_.p_ *(params_.p_ + 2), MemcpyHostToDevice);

    idx = 0;
    for(int ii=0; ii< 2 * params_.rep_up_freq_; ++ii)
    {
        len = params_.rep_up_freq_ + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
            buffer[idx++] = (len-jj)<=(params_.rep_up_freq_-params_.rep_filter_freq_) ? 0 : 1;
    }
    device_.Memcpy(buffer, alpha_q,params_.rep_up_freq_ *(params_.rep_up_freq_ + 2), MemcpyHostToDevice);

    free(buffer);
    SetX(x_in);
}

template<typename T>
Surface<T>::~Surface()
{
    device_.Free(shc);
    device_.Free(work_arr);
    device_.Free(alpha_p);
    device_.Free(alpha_q);
    device_.Free(tension_);
    device_.Free(all_rot_mats_);
    device_.Free(rot_mat);
    device_.Free(sing_quad_weights_);
    device_.Free(vel);
    device_.Free(quad_weights_);
}

template<typename T>
void Surface<T>::Resize(int n_surfs_in)
{
    x_.Resize(n_surfs_in);
    normal_.Resize(n_surfs_in); 
    h_.Resize(n_surfs_in); 
    w_.Resize(n_surfs_in);
    k_.Resize(n_surfs_in); 
    cu_.Resize(n_surfs_in); 
    cv_.Resize(n_surfs_in);
    bending_force_.Resize(n_surfs_in);
    tensile_force_.Resize(n_surfs_in);
    S1.Resize(n_surfs_in);
    S2.Resize(n_surfs_in);
    S3.Resize(n_surfs_in);
    S4.Resize(n_surfs_in); 
    S5.Resize(n_surfs_in); 
    S6.Resize(n_surfs_in);
    V1.Resize(n_surfs_in);
    V2.Resize(n_surfs_in);
    V10.Resize(n_surfs_in); 
    V11.Resize(n_surfs_in);
    V12.Resize(n_surfs_in); 
    V13.Resize(n_surfs_in);
    S10.Resize(n_surfs_in);
    w_sph_(n_surfs_in);

}
template <typename T> 
void Surface<T>::SetX(const Vectors<T> &x_in)
{
    x_.SetData(x_in.data_);
    UpdateAll();
}

template <typename T> 
void Surface<T>::UpdateFirstForms()
{
    // Spherical harmonic coefficient
    device_.ShAna(x_.data_, work_arr, params_.p_, 3*params_.n_surfs_, shc);

    // First derivatives
    device_.ShSynDu(shc, work_arr, params_.p_, 3*params_.n_surfs_, V1.data_);//V1 is holding xu
    device_.ShSynDv(shc, work_arr, params_.p_, 3*params_.n_surfs_, V2.data_);//V2 is holding xv

    // First fundamental coefficients
    DotProduct(V1,V1, S1);
    DotProduct(V1,V2, S2);
    DotProduct(V2,V2, S3);
    CrossProduct(V1,V2,normal_);

    // Area element
    DotProduct(normal_, normal_, w_);//w = W^2 for now
    
    // Dividing EFG by W^2
    xyInv(S1, w_, S1);
    xyInv(S2, w_, S2);
    xyInv(S3, w_, S3);
    w_.Sqrt();
    uyInv(normal_,w_,normal_);
}

template <typename T> 
void Surface<T>::UpdateAll()
{
    UpdateFirstForms();///@todo This may lead to a bug. Since there is no guarantee that E,S2,S3 hold their values.

    //Div and Grad coefficients
    xvpb(S2, V2, (T) 0.0, cu_);//F*xv
    axpb((T) -1.0,cu_, (T) 0.0, cu_);
    xvpw(S3, V1, cu_, cu_);
     
    xvpb(S2, V1, (T) 0.0, cv_);//F*xu
    axpb((T) -1.0,cv_, (T) 0.0, cv_);
    xvpw(S1, V2, cv_, cv_);

    // Second derivatives
    device_.ShSynDuu(shc, work_arr, params_.p_, 3*params_.n_surfs_, V1.data_);//V1 is holding xuu
    DotProduct(V1, normal_, S4);
    
    device_.ShSynDuv(shc, work_arr, params_.p_, 3*params_.n_surfs_, V1.data_);//V1 is holding xuv
    DotProduct(V1, normal_, S5);

    device_.ShSynDvv(shc, work_arr, params_.p_, 3*params_.n_surfs_, V1.data_);//V1 is holding xvv
    DotProduct(V1, normal_, S6);
    
    // Gaussian curvature
    xy(S4,S6,k_);
    xy(S5,S5,h_);
    axpy((T)-1.0,h_,k_,k_);
    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);
    
    device_.Filter(params_.p_, params_.n_surfs_, k_.data_, alpha_p, work_arr, shc, k_.data_);
    
    // Mean curvature
    xy(S1,S6,h_);
    axpb((T) .5,h_, (T) 0.0,h_);
    
    xy(S2,S5,S6);
    axpy((T)-1.0, S6, h_, h_);
    
    xy(S3,S4,S6);
    axpb((T) .5,S6, (T) 0.0, S6);
    
    axpy((T)1.0 ,h_,S6, h_);

    device_.Filter(params_.p_, params_.n_surfs_, h_.data_, alpha_p, work_arr, shc, h_.data_);

    //Bending force
    SurfGrad(h_, bending_force_);
    SurfDiv(bending_force_, S1);
    axpb((T) -1.0,S1, (T) 0.0, S1);
    xy(h_,h_,S2);
    axpy((T) -1.0, S2, k_, S2);
    xy(S2, h_, S2);
    axpy((T) 2.0, S2, S1, S1);
    xvpb(S1, normal_, (T) 0.0, bending_force_);
    axpb(params_.kappa_,bending_force_, (T) 0.0, bending_force_);
    
    //Tensile force
    xvpb(h_,normal_,(T) 0.0, tensile_force_);
}

template <typename T>
void Surface<T>::StokesMatVec(const Vectors<T> &density_in, Vectors<T> &velocity_out)
{
#ifdef PROFILING
    double ss = get_seconds();
#endif
        
    int np = 2 * params_.p_ * (params_.p_ + 1);
    int rot_chunck = 2 * params_.p_ * np;
    T alpha(1.0), beta(0.0);
    int trg_idx(0);
    
    DataIO<T> myIO;
    char fname[300];
    
    xyInv(w_, w_sph_, S1);
   
    //@todo #pragma omp parallel for 
    int nvX3= 3*params_.n_surfs_;
    for(int ii=0;ii <= params_.p_; ++ii)
    {
        for(int jj=0;jj < 2 * params_.p_; ++jj)
        {
            device_.CircShift( all_rot_mats_ + ii * np * np, params_.p_ + 1, rot_chunck, jj * np, rot_mat);
            
            device_.gemm("N", "N", &np, &nvX3, &np, 
			 &alpha,rot_mat, &np, x_.data_, &np, &beta, V1.data_, &np);
            
            device_.gemm("N", "N", &np, &nvX3, &np,
			 &alpha,rot_mat,&np,density_in.data_,&np,&beta, V2.data_, &np);
            
            device_.gemm("N", "N", &np,&params_.n_surfs_, &np,
			 &alpha, rot_mat, &np, S1.data_, &np, &beta, S2.data_, &np);
    
            xy(S2, w_sph_, S2);
            xvpb(S2, V2, (T) 0.0, V2);
            device_.DirectStokes(np, params_.n_surfs_, ii*2*params_.p_ + jj, ii*2*params_.p_ + jj + 1, 
                sing_quad_weights_, x_.data_, V1.data_, V2.data_, velocity_out.data_);
            
        }
    }
#ifdef PROFILING
    ss = get_seconds()-ss ;
    cout<<"StokesMatVec (sec) : "<<ss<<endl;
#endif

}

template <typename T> 
void Surface<T>::SurfGrad(const Scalars<T> &f_in, Vectors<T> &grad_f_out)
{
    device_.FirstDerivatives(f_in.data_, work_arr, params_.p_, params_.n_surfs_, shc, S5.data_, S6.data_);
    
    xvpb(S5,cu_, (T) 0.0, grad_f_out);
    xvpw(S6,cv_, grad_f_out, grad_f_out);

    device_.Filter(params_.p_, 3*params_.n_surfs_, grad_f_out.data_, alpha_p, work_arr, shc, grad_f_out.data_);
}

///@todo this can be done in place so that we'll only need one extra work space.
template <typename T> 
void Surface<T>::SurfDiv(const Vectors<T> &f_in, Scalars<T> &div_f_out) 
{
    device_.FirstDerivatives(f_in.data_, work_arr, params_.p_, 3*params_.n_surfs_, shc, V1.data_, V2.data_);
    
    DotProduct(V1,cu_, S6);
    DotProduct(V2,cv_, div_f_out);
    axpy((T)1.0,div_f_out,S6, div_f_out);

    device_.Filter(params_.p_, params_.n_surfs_, div_f_out.data_, alpha_p, work_arr, shc, div_f_out.data_);
}

template <typename T> 
void Surface<T>::Reparam()
{
    int iter = 0;
    T vel;

    device_.ShAna(x_.data_, work_arr, params_.p_, 3*params_.n_surfs_, V10.data_);
    device_.Resample(params_.p_, 3*params_.n_surfs_, params_.rep_up_freq_, V10.data_, shc);
    device_.ShSyn(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V10.data_);

    //upsample
    while(iter++ < params_.rep_iter_max_)
    {
        UpdateNormal();///@todo shc may have changed
        device_.Filter(params_.rep_up_freq_, 3*params_.n_surfs_, V10.data_, alpha_q, work_arr, shc, V11.data_);

        axpy((T) -1.0, V10, V11, V11);
        DotProduct(V13, V11, S10);
        axpb((T) -1.0,S10, (T) 0.0, S10);
        xvpw(S10, V13, V11, V11);//The correction velocity
        axpy(params_.rep_ts_, V11, V10, V10);

#ifndef NDEBUG
        DotProduct(V11,V11,S10);//@todo seems extra because we work with singles now
        vel = S10.Max();
        cout<<" Iteration"<<iter<<", reparam max vel: "<<vel<<endl;
#endif
    }
    device_.ShAna(V10.data_, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V11.data_);
    device_.Resample(params_.rep_up_freq_, 3*params_.n_surfs_, params_.p_, V11.data_, shc);
    device_.ShSyn(shc, work_arr, params_.p_, 3*params_.n_surfs_, x_.data_);
}

template<typename T>
void Surface<T>::UpdateNormal()
{
    device_.ShAna(V10.data_, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, shc);
    device_.ShSynDu(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V11.data_);
    device_.ShSynDv(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V12.data_);

    CrossProduct(V11,V12,V13);
    DotProduct(V13, V13, S10);
    S10.Sqrt();
    uyInv(V13,S10,V13);
}

template<typename T>
void Surface<T>::Area()
{
    device_.Reduce(NULL, w_.data_, quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
    for(int ii=0;ii<params_.n_surfs_;++ii)
        cout<<work_arr[ii]<<endl;
}

template<typename T>
void Surface<T>::Volume()
{
    DotProduct(x_,normal_,S1);
    axpb((T) 1/3,S1,(T) 0.0, S1);
    device_.Reduce(S1.data_, w_.data_, quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
    for(int ii=0;ii<params_.n_surfs_;++ii)
        cout<<work_arr[ii]<<endl;
}

template<typename T>
void Surface<T>::GetTension(const Vectors<T> &v_in, const Vectors<T> &v_ten_in, T *tension_out)
{
    SurfDiv(v_in, S1);
    SurfDiv(v_ten_in, S2);
    
    device_.Reduce(S1.data_, w_.data_, quad_weights_, S1.GetFunLength(), params_.n_surfs_, tension_out);
    device_.Reduce(S2.data_, w_.data_, quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
    device_.axpb((T) -1.0, tension_out, (T) 0.0, 1, params_.n_surfs_, tension_out);
    device_.xyInv(tension_out, work_arr, 1, params_.n_surfs_, tension_out);
}
