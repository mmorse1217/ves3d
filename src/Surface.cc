/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */

template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::Surface(SurfaceParams<value_type> &params_in, 
    const OperatorsMats<value_type> &mats) :
    params_(params_in),
    x_(params_.p_,params_.n_surfs_),
    normal_( params_.p_,params_.n_surfs_), 
    h_(params_.p_,params_.n_surfs_), 
    w_(params_.p_,params_.n_surfs_),
    k_(params_.p_,params_.n_surfs_), 
    cu_(params_.p_,params_.n_surfs_), 
    cv_(params_.p_,params_.n_surfs_),
    //bending_force_(params_.p_,params_.n_surfs_),
    //tension_( params_.p_, params_.n_surfs_, make_pair(1,1)),
    //tensile_force_(params_.p_,params_.n_surfs_),
    capacity_(params_.n_surfs_),
    //GB: the code below is ugly, must be fixed after sc10
    S1(params_.p_,params_.n_surfs_),
    S2(params_.p_,params_.n_surfs_),
    S3(params_.p_,params_.n_surfs_),
    S4(params_.p_,params_.n_surfs_), 
    S5(params_.p_,params_.n_surfs_), 
    S6(params_.p_,params_.n_surfs_),
    V1(params_.p_,params_.n_surfs_),
    V2(params_.p_,params_.n_surfs_),
//     V10(params_.rep_up_freq_,params_.n_surfs_), 
//     V11(params_.rep_up_freq_,params_.n_surfs_),
//     V12(params_.rep_up_freq_,params_.n_surfs_),
//     V13(params_.rep_up_freq_,params_.n_surfs_),
//     S10(params_.rep_up_freq_,params_.n_surfs_),
    //w_sph_(params_.p_,params_.n_surfs_),
    //max_init_area_(-1),
    //StokesMatVec_time_(0),
    sht_(&x_.GetDevice(), params_.p_)
{
    shc = (value_type*) x_.GetDevice().Malloc(6  * params_.rep_up_freq_ * 
        (params_.rep_up_freq_ + 2) * params_.n_surfs_ * sizeof(value_type));
    
    work_arr = (value_type*) x_.GetDevice().Malloc(12 * params_.rep_up_freq_ * 
        (params_.rep_up_freq_ + 1) * params_.n_surfs_ * sizeof(value_type));
    
    //alpha_p  = (T*) device_->Malloc(params_.p_ *(params_.p_ + 2) * sizeof(T));
    //alpha_q  = (T*) device_->Malloc(params_.rep_up_freq_ *(params_.rep_up_freq_ + 2) * sizeof(T));
    //tension_ = (T*) device_->Calloc(params_.n_surfs_, sizeof(T));

    //int np = 2 * params_.p_ * (params_.p_ + 1);

    //quad_weights_ = mats.quad_weights_;
    //all_rot_mats_ = mats.all_rot_mats_;
    //sing_quad_weights_ = mats.sing_quad_weights_;
    
    //device_->Memcpy(w_sph_.begin(), mats.w_sph_, np * sizeof(T), MemcpyDeviceToDevice);
    //for(int ii=1;ii<params_.n_surfs_;++ii)
    //    device_->Memcpy(w_sph_.begin() + ii*np, w_sph_.begin(), np * sizeof(T), MemcpyDeviceToDevice);

    //rot_mat = (T*) device_->Malloc(np * np * sizeof(T));
    //T *buffer = (T*) malloc(np * np * (params_.p_ + 1) * sizeof(T));

    //int idx = 0, len;
    //for(int ii=0; ii< 2 * params_.p_; ++ii)
    //{
    //     len = params_.p_ + 1 - (ii+1)/2;
//         for(int jj=0; jj < len; ++jj)
//             buffer[idx++] = (len-jj)<=(params_.p_-params_.filter_freq_) ? 0 : 1;
//     }
//     device_->Memcpy(alpha_p, buffer, params_.p_ *(params_.p_ + 2) * sizeof(T), MemcpyHostToDevice);

//     idx = 0;
//     for(int ii=0; ii< 2 * params_.rep_up_freq_; ++ii)
//     {
//         len = params_.rep_up_freq_ + 1 - (ii+1)/2;
//         for(int jj=0; jj < len; ++jj)
//             buffer[idx++] = (len-jj)<=(params_.rep_up_freq_-params_.rep_filter_freq_) ? 0 : 1;
//     }

//     device_->Memcpy(alpha_q, buffer, params_.rep_up_freq_ *(params_.rep_up_freq_ + 2) * sizeof(T), MemcpyHostToDevice);

//     free(buffer);
}

template<typename ScalarContainer, typename VectorContainer> 
Surface<ScalarContainer, VectorContainer>::~Surface()
{
    x_.GetDevice().Free(shc);
    x_.GetDevice().Free(work_arr);
    //device_->Free(alpha_p);
    //device_->Free(alpha_q);
    //device_->Free(tension_);
    //device_->Free(rot_mat);
}

template<typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::Resize(int n_surfs_in)
{
    x_.Resize(n_surfs_in);
    normal_.Resize(n_surfs_in); 
    cu_.Resize(n_surfs_in); 
    cv_.Resize(n_surfs_in);
    w_.Resize(n_surfs_in);
    h_.Resize(n_surfs_in); 
    k_.Resize(n_surfs_in); 
    //bending_force_.Resize(n_surfs_in);
    //tensile_force_.Resize(n_surfs_in);
    S1.Resize(n_surfs_in);
    S2.Resize(n_surfs_in);
    S3.Resize(n_surfs_in);
    S4.Resize(n_surfs_in); 
    S5.Resize(n_surfs_in); 
    S6.Resize(n_surfs_in);
    V1.Resize(n_surfs_in);
    V2.Resize(n_surfs_in);
    // S10.Resize(n_surfs_in);
//     V10.Resize(n_surfs_in); 
//     V11.Resize(n_surfs_in);
//     V12.Resize(n_surfs_in); 
//     V13.Resize(n_surfs_in);
//     tension_.Resize(n_surfs_in, params_.p_, make_pair(1,1));
    
#ifndef NDEBUG
    cout<<"Resizing form "<<params_.n_surfs_<<" to "<<n_surfs_in<<endl;
#endif

    if(n_surfs_in > capacity_)
    {
#ifndef NDEBUG
        cout<<"  . The new size is larger than the current allocated memory, allocating new memory. "<<endl;
        cout<<"  . New size "<<n_surfs_in<<endl;
#endif
        this->capacity_ = n_surfs_in;
        x_.GetDevice().Free(shc);
        shc = (value_type*) x_.GetDevice().Malloc(6  * params_.rep_up_freq_ * (params_.rep_up_freq_ + 1) * capacity_ * sizeof(value_type));
        
        x_.GetDevice().Free(work_arr);
        work_arr = (value_type*) x_.GetDevice().Malloc(12 * params_.rep_up_freq_ * (params_.rep_up_freq_ + 1) * capacity_ * sizeof(value_type));
        
//         T *tension_old(this->tension_);
//         tension_ = (T*) x_.GetDevice().Calloc(capacity_, sizeof(T));
//         if(tension_old != 0)
//             x_.GetDevice().Memcpy(tension_, tension_old, params_.n_surfs_ * sizeof(T), MemcpyDeviceToDevice);
//         x_.GetDevice().Free(tension_old);
        
        //w_sph_.Resize(capacity_);
        //int np = 2 * params_.p_ * (params_.p_ + 1);
        //for(int ii=params_.n_surfs_;ii<capacity_;++ii)
        //   x_.GetDevice().Memcpy(w_sph_.begin() + ii*np, w_sph_.begin(), np * sizeof(T), MemcpyDeviceToDevice);
    }
    this->params_.n_surfs_ = n_surfs_in;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::UpdateFirstForms()
{
    // Spherical harmonic coefficient
    sht_.forward(x_.begin(), work_arr, 3*params_.n_surfs_, shc);

    // First derivatives
    sht_.backward_du(shc, work_arr, 3*params_.n_surfs_, V1.begin());
    sht_.backward_dv(shc, work_arr, 3*params_.n_surfs_, V2.begin());

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
    Sqrt(w_, w_);
    uyInv(normal_,w_,normal_);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::UpdateAll()
{
    UpdateFirstForms();///@todo This may lead to a bug. Since there is no guarantee that S1,S2,S3 hold their values.

    //Div and Grad coefficients
    xv(S2, V2, cu_);//F*xv
    axpy((value_type) -1.0,cu_, cu_);
    xvpw(S3, V1, cu_, cu_);
     
    xv(S2, V1, cv_);//F*xu
    axpy((value_type) -1.0,cv_, cv_);
    xvpw(S1, V2, cv_, cv_);

    // Second derivatives
    sht_.backward_d2u(shc, work_arr, 3*params_.n_surfs_, V1.begin());
    DotProduct(V1, normal_, S4);
    
    sht_.backward_duv(shc, work_arr, 3*params_.n_surfs_, V1.begin());
    DotProduct(V1, normal_, S5);

    sht_.backward_d2v(shc, work_arr, 3*params_.n_surfs_, V1.begin());
    DotProduct(V1, normal_, S6);
    
    // Gaussian curvature
    xy(S4,S6,k_);
    xy(S5,S5,h_);
    axpy((value_type)-1.0,h_,k_,k_);
    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);
    
    //x_.GetDevice().Filter(params_.p_, params_.n_surfs_, k_.begin(), alpha_p, work_arr, shc, k_.begin());
    
    // Mean curvature
    xy(S1,S6,h_);
    axpy((value_type) .5,h_,h_);
    
    xy(S2,S5,S6);
    axpy((value_type)-1.0, S6, h_, h_);
    
    xy(S3,S4,S6);
    axpy((value_type) .5,S6, S6);
    
    axpy((value_type)1.0 ,h_,S6, h_);

    //x_.GetDevice().Filter(params_.p_, params_.n_surfs_, h_.begin(), alpha_p, work_arr, shc, h_.begin());

//     //Bending force
//     Grad(h_, bending_force_);
//     Div(bending_force_, S1);
//     axpy((value_type) -1.0,S1, S1);
//     xy(h_,h_,S2);
//     axpy((value_type) -1.0, S2, k_, S2);
//     xy(S2, h_, S2);
//     axpy((value_type) 2.0, S2, S1, S1);
//     xv(S1, normal_, bending_force_);
//     axpy(params_.kappa_,bending_force_, bending_force_);
    
//     //Tensile force
//     xv(h_,normal_, tensile_force_);

    //if(max_init_area_<0)
    //    max_init_area_ = this->Area();
}

// template <typename ScalarContainer, typename VectorContainer> 
// void Surface<ScalarContainer, VectorContainer>::StokesMatVec(const VectorContainer &density_in, VectorContainer &velocity_out)
// {
//     int np = 2 * params_.p_ * (params_.p_ + 1);
//     int rot_chunck = 2 * params_.p_ * np;
//     T alpha(1.0), beta(0.0);
//     int trg_idx(0);
    
//     xyInv(w_, w_sph_, S1);
//     int nvX3= 3*params_.n_surfs_;   

//     for(int ii=0;ii <= params_.p_; ++ii)
//     {
//         for(int jj=0;jj < 2 * params_.p_; ++jj)
//         {
//             x_.GetDevice().CircShift( all_rot_mats_ + ii * np * np, params_.p_ + 1, rot_chunck, jj * np, rot_mat);
            
//             x_.GetDevice().gemm("N", "N", &np, &nvX3, &np, 
//                 &alpha,rot_mat, &np, x_.begin(), &np, &beta, V1.begin(), &np);
    
//             x_.GetDevice().gemm("N", "N", &np, &nvX3, &np,
//                 &alpha,rot_mat,&np,density_in.begin(),&np,&beta, V2.begin(), &np);
            
//             x_.GetDevice().gemm("N", "N", &np,&params_.n_surfs_, &np,
//                 &alpha, rot_mat, &np, S1.begin(), &np, &beta, S2.begin(), &np);
    
//             xy(S2, w_sph_, S2);
//             xvpb(S2, V2, (value_type) 0.0, V2);

//             x_.GetDevice().DirectStokes(np, params_.n_surfs_, ii*2*params_.p_ + jj, ii*2*params_.p_ + jj + 1, 
//                 sing_quad_weights_, x_.begin(), V1.begin(), V2.begin(), velocity_out.begin());
            
//         }
//     }

// }

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::Grad(const ScalarContainer &f_in, VectorContainer &grad_f_out) const
{
    sht_.forward( f_in.begin(), work_arr, params_.n_surfs_, shc);
    sht_.backward_du( shc, work_arr, params_.n_surfs_, S5.begin());
    sht_.backward_dv( shc, work_arr, params_.n_surfs_, S6.begin());

    xv(S5,cu_, grad_f_out);
    xvpw(S6,cv_, grad_f_out, grad_f_out);

    //x_.GetDevice().Filter(params_.p_, 3*params_.n_surfs_, grad_f_out.begin(), alpha_p, work_arr, shc, grad_f_out.begin());
}

///@todo this can be done in place so that we'll only need one extra work space.
template <typename ScalarContainer, typename VectorContainer > 
void Surface<ScalarContainer, VectorContainer>::Div(const VectorContainer &f_in, ScalarContainer &div_f_out) const
{
    sht_.forward( f_in.begin(), work_arr, 3*params_.n_surfs_, shc);
    sht_.backward_du( shc, work_arr, 3*params_.n_surfs_, V1.begin());
    sht_.backward_dv( shc, work_arr, 3*params_.n_surfs_, V2.begin());
    
    DotProduct(V1,cu_, S6);
    DotProduct(V2,cv_, div_f_out);
    axpy((value_type)1.0,div_f_out,S6, div_f_out);

    //x_.GetDevice().Filter(params_.p_, params_.n_surfs_, div_f_out.begin(), alpha_p, work_arr, shc, div_f_out.begin());
}

// template <typename ScalarContainer, typename VectorContainer > 
// void Surface<ScalarContainer, VectorContainer>::Reparam()
// {
//     int iter = 0;
//     T vel = 2*params_.rep_max_vel_;

//     //up-sample
//     x_.GetDevice().InterpSh(params_.p_, 3*params_.n_surfs_, x_.begin(), work_arr, shc, params_.rep_up_freq_, V10.begin());

//     while(iter++ < params_.rep_iter_max_ && vel > params_.rep_max_vel_ * params_.rep_max_vel_)
//     {
//         UpdateNormal();///@todo shc may have changed
//         x_.GetDevice().Filter(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), alpha_q, work_arr, shc, V11.begin());

//         axpy((T) -1.0, V10, V11, V11);
//         DotProduct(V13, V11, S10);
//         axpb((T) -1.0,S10, (T) 0.0, S10);
//         xvpw(S10, V13, V11, V11);//The correction velocity
//         axpy(params_.rep_ts_, V11, V10, V10);

//         DotProduct(V11,V11,S10);
//         ///@todo The error should be relative
//         vel = S10.Max();
// #ifndef NDEBUG
//         cout<<" Reparam iteration "<<iter<<", max vel: "<<vel<<endl;
// #endif
//     }
    
//     ///inteprsh is buggy
//     x_.GetDevice().InterpSh(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), work_arr, shc, params_.p_, V11.begin());
//     x_.GetDevice().Memcpy(x_.begin(), V11.begin(), x_.GetDataLength(), MemcpyDeviceToDevice);
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::UpdateNormal()
// {
//     x_.GetDevice().ShAna(V10.begin(), work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, shc);
//     x_.GetDevice().ShSynDu(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V11.begin());
//     x_.GetDevice().ShSynDv(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V12.begin());

//     CrossProduct(V11,V12,V13);
//     DotProduct(V13, V13, S10);
//     S10.Sqrt();
//     uyInv(V13,S10,V13);
// }

// template< typename ScalarContainer, typename VectorContainer >
// T Surface<ScalarContainer, VectorContainer>::Area()
// {
//     //Reduce(w_.begin(), quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
//     //return Max(work_arr, params_.n_surfs_);
//     return 0;
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::Volume()
// {
//     DotProduct(x_,normal_,S1);
//     axpy((T) 1/3,S1, S1);
//     Reduce(S1.begin(), w_.begin(), quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
//     //     for(int ii=0;ii<params_.n_surfs_;++ii)
//     //         cout<<work_arr[ii]<<endl;
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::GetTension(const VectorContainer &v_in, const VectorContainer &v_ten_in, T *tension_out)
// {
//     SurfDiv(v_in, S1);
//     SurfDiv(v_ten_in, S2);
    
//     x_.GetDevice().Reduce(S1.begin(), w_.begin(), quad_weights_, S1.GetFunLength(), params_.n_surfs_, tension_out);
//     x_.GetDevice().Reduce(S2.begin(), w_.begin(), quad_weights_, S1.GetFunLength(), params_.n_surfs_, work_arr);
//     x_.GetDevice().axpb((T) -1.0, tension_out, (T) 0.0, 1, params_.n_surfs_, tension_out);
//     x_.GetDevice().xyInv(tension_out, work_arr, 1, params_.n_surfs_, tension_out);
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::Populate(const T *centers)
// {
//     int length = this->x_.GetFunLength();
//     for(int ii=1;ii<params_.n_surfs_;ii++)
//     {
//         int idx = 3 * ii * length;
//         x_.GetDevice().Memcpy(x_.begin() + idx, x_.begin(), 3 * length, MemcpyDeviceToDevice);
//         x_.GetDevice().axpb((T) 1.0, x_.begin() + idx                  , centers[3*ii  ], length, 1, x_.begin() + idx                  );
//         x_.GetDevice().axpb((T) 1.0, x_.begin() + idx + length         , centers[3*ii+1], length, 1, x_.begin() + idx + length         );
//         x_.GetDevice().axpb((T) 1.0, x_.begin() + idx + length + length, centers[3*ii+2], length, 1, x_.begin() + idx + length + length);
//     }

//     //treating the first surface
//     x_.GetDevice().axpb((T) 1.0, x_.begin()                  , centers[0], length, 1, x_.begin()                  );
//     x_.GetDevice().axpb((T) 1.0, x_.begin() + length         , centers[1], length, 1, x_.begin() + length         );
//     x_.GetDevice().axpb((T) 1.0, x_.begin() + length + length, centers[2], length, 1, x_.begin() + length + length);
// }

// template<typename ScalarContainer, typename VectorContainer >
// T* Surface<ScalarContainer, VectorContainer>::GetCenters(T* cnts)
// {
//     UpdateFirstForms();
    
//     DotProduct(x_, x_, S1);
//     axpb((T) 1/2, S1,(T) 0.0, S1);
//     xvpb(S1, normal_, (T) 0.0, V1);
    
//     size_t idx = 0;
//     size_t len = w_.GetDataLength();
//     size_t sc_len = w_.GetFunLength();
//     size_t nv = params_.n_surfs_;

//     x_.GetDevice().Memcpy(work_arr            , w_.begin(), len, MemcpyDeviceToDevice);
//     x_.GetDevice().Memcpy(work_arr + len      , w_.begin(), len, MemcpyDeviceToDevice);
//     x_.GetDevice().Memcpy(work_arr + len + len, w_.begin(), len, MemcpyDeviceToDevice);

//     x_.GetDevice().ShufflePoints(work_arr                   , AxisMajor , len   , 1 , work_arr + len + len + len);
//     x_.GetDevice().ShufflePoints(work_arr + len + len + len , PointMajor, sc_len, nv, work_arr);

//     x_.GetDevice().Reduce(NULL, work_arr, quad_weights_, sc_len, 3*nv, cnts);
//     x_.GetDevice().Reduce(V1.begin(), work_arr, quad_weights_, sc_len, 3*nv, cnts);

//     DotProduct(x_,normal_,S1);
//     axpb((T) 1/3,S1,(T) 0.0, S1);
//     x_.GetDevice().Reduce(S1.begin(), w_.begin(), quad_weights_, sc_len, nv, work_arr);
    
//     x_.GetDevice().Memcpy(work_arr + nv     , work_arr, nv, MemcpyDeviceToDevice); 
//     x_.GetDevice().Memcpy(work_arr + nv + nv, work_arr, nv, MemcpyDeviceTDevice); 
//     x_.GetDevice().ShufflePoints(work_arr, AxisMajor, nv, 1, work_arr + 3*nv);
//     x_.GetDevice().xyInv(cnts, work_arr + 3*nv, 3, nv, cnts);

//     return cnts;
// }
