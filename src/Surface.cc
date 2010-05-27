/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */
template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::Surface(const Vec& x_in) :
    sht_(&x_.getDevice(), x_in.getShOrder()), ///@todo fix this
    position_has_changed_outside_(true),
    first_forms_are_stale_(true),
    second_forms_are_stale_(true)
{
    setPosition(x_in);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::setPosition(const Vec& x_in)
{
    position_has_changed_outside_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    x_.replicate(x_in);
    axpy(1.0, x_in, x_);
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer& Surface<ScalarContainer, VectorContainer>::getPosition()
{
    position_has_changed_outside_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, 
                               VectorContainer>::getPosition() const
{
    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, 
                               VectorContainer>::getNormal() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(normal_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getAreaElement() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(w_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getMeanCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(h_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, 
                               VectorContainer>::getGaussianCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(k_);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateFirstForms() const
{
    if(position_has_changed_outside_)
        checkContainers();
        
    size_t nv3 = x_.getTheDim() * x_.getNumSubs();
    // Spherical harmonic coefficient
    sht_.forward(x_.begin(), work_arr.begin(), nv3, shc.begin());
    
    // First derivatives
    sht_.backward_du(shc.begin(), work_arr.begin(), nv3, V1.begin());
    sht_.backward_dv(shc.begin(), work_arr.begin(), nv3, V2.begin());
    
    // First fundamental coefficients
    GeometricDot(V1,V1, S1);
    GeometricDot(V1,V2, S2);
    GeometricDot(V2,V2, S3);
    GeometricCross(V1,V2,normal_);
    
    // Area element
    GeometricDot(normal_, normal_, w_);//w = W^2 for now
    
    // Dividing EFG by W^2
    xyInv(S1, w_, S1);
    xyInv(S2, w_, S2);
    xyInv(S3, w_, S3);
    Sqrt(w_, w_);
    uyInv(normal_,w_,normal_);

    first_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::checkContainers() const
{
    normal_.replicate(x_);
    w_.replicate(x_);    
    h_.replicate(x_);    
    k_.replicate(x_);    
    cu_.replicate(x_);    
    cv_.replicate(x_);    
    
    ///@todo change these
    shc.replicate(x_);
    work_arr.resize(2 * x_.getNumSubs(), x_.getShOrder(), x_.getGridDim());
    S1.replicate(x_);    
    S2.replicate(x_);    
    S3.replicate(x_);    
    S4.replicate(x_);    
    S5.replicate(x_);    
    S6.replicate(x_);    
    V1.replicate(x_);    
    V2.replicate(x_);

}

// {
//     shc = (value_type*) x_.getDevice().Malloc(6  * params_.rep_up_freq_ * 
//         (params_.rep_up_freq_ + 2) * params_.n_surfs_ * sizeof(value_type));
    
//     work_arr = (value_type*) x_.getDevice().Malloc(12 * params_.rep_up_freq_ * 
//         (params_.rep_up_freq_ + 1) * params_.n_surfs_ * sizeof(value_type));   

//template<typename ScalarContainer, typename VectorContainer> 
// Surface<ScalarContainer, VectorContainer>::~Surface()
// {
//     x_.getDevice().Free(shc);
//     x_.getDevice().Free(work_arr);
// }

// template<typename ScalarContainer, typename VectorContainer> 
// void Surface<ScalarContainer, VectorContainer>::Resize(int n_surfs_in)
// {
//     x_.resize(n_surfs_in);
//     normal_.resize(n_surfs_in); 
//     cu_.resize(n_surfs_in); 
//     cv_.resize(n_surfs_in);
//     w_.resize(n_surfs_in);
//     h_.resize(n_surfs_in); 
//     k_.resize(n_surfs_in); 
//     S1.resize(n_surfs_in);
//     S2.resize(n_surfs_in);
//     S3.resize(n_surfs_in);
//     S4.resize(n_surfs_in); 
//     S5.resize(n_surfs_in); 
//     S6.resize(n_surfs_in);
//     V1.resize(n_surfs_in);
//     V2.resize(n_surfs_in);
    
// #ifndef NDEBUG
//     cout<<"Resizing form "<<params_.n_surfs_<<" to "<<n_surfs_in<<endl;
// #endif

//     if(n_surfs_in > capacity_)
//     {
// #ifndef NDEBUG
//         cout<<"  . The new size is larger than the current allocated memory, allocating new memory. "<<endl;
//         cout<<"  . New size "<<n_surfs_in<<endl;
// #endif
//         this->capacity_ = n_surfs_in;
//         x_.getDevice().Free(shc);
//         shc = (value_type*) x_.getDevice().Malloc(6  * params_.rep_up_freq_ * (params_.rep_up_freq_ + 1) * capacity_ * sizeof(value_type));
        
//         x_.getDevice().Free(work_arr);
//         work_arr = (value_type*) x_.getDevice().Malloc(12 * params_.rep_up_freq_ * (params_.rep_up_freq_ + 1) * capacity_ * sizeof(value_type));

//     }
//     this->params_.n_surfs_ = n_surfs_in;
// }

// template <typename ScalarContainer, typename VectorContainer>  
// void Surface<ScalarContainer, VectorContainer>::UpdateAll()
// {
//     ///@todo This may lead to a bug. Since there is no guarantee that
//     ///S1,S2,S3 hold their values.
//     UpdateFirstForms();

//     //Div and Grad coefficients
//     xv(S2, V2, cu_);//F*xv
//     axpy((value_type) -1.0,cu_, cu_);
//     xvpw(S3, V1, cu_, cu_);
     
//     xv(S2, V1, cv_);//F*xu
//     axpy((value_type) -1.0,cv_, cv_);
//     xvpw(S1, V2, cv_, cv_);

//     // Second derivatives
//     sht_.backward_d2u(shc, work_arr, 3*params_.n_surfs_, V1.begin());
//     GeometricDot(V1, normal_, S4);
    
//     sht_.backward_duv(shc, work_arr, 3*params_.n_surfs_, V1.begin());
//     GeometricDot(V1, normal_, S5);

//     sht_.backward_d2v(shc, work_arr, 3*params_.n_surfs_, V1.begin());
//     GeometricDot(V1, normal_, S6);
    
//     // Gaussian curvature
//     xy(S4,S6,k_);
//     xy(S5,S5,h_);
//     axpy((value_type)-1.0,h_,k_,k_);
//     xyInv(k_,w_,k_);
//     xyInv(k_,w_,k_);
    
//     sht_.Filter(k_.begin(),work_arr, params_.n_surfs_, shc, k_.begin());
    
//     // Mean curvature
//     xy(S1,S6,h_);
//     axpy((value_type) .5,h_,h_);
    
//     xy(S2,S5,S6);
//     axpy((value_type)-1.0, S6, h_, h_);
    
//     xy(S3,S4,S6);
//     axpy((value_type) .5,S6, S6);
    
//     axpy((value_type)1.0 ,h_,S6, h_);

//     sht_.Filter(h_.begin(),work_arr, params_.n_surfs_, shc, h_.begin());
// }

// template <typename ScalarContainer, typename VectorContainer>  
// void Surface<ScalarContainer, VectorContainer>::Grad(const ScalarContainer &f_in, VectorContainer &grad_f_out) const
// {
//     sht_.forward( f_in.begin(), work_arr, params_.n_surfs_, shc);
//     sht_.backward_du( shc, work_arr, params_.n_surfs_, S5.begin());
//     sht_.backward_dv( shc, work_arr, params_.n_surfs_, S6.begin());

//     xv(S5,cu_, grad_f_out);
//     xvpw(S6,cv_, grad_f_out, grad_f_out);
//     sht_.Filter(grad_f_out.begin(), work_arr, 3*params_.n_surfs_, shc, grad_f_out.begin());
// }

// ///@todo this can be done in place so that we'll only need one extra work space.
// template <typename ScalarContainer, typename VectorContainer > 
// void Surface<ScalarContainer, VectorContainer>::Div(const VectorContainer &f_in, ScalarContainer &div_f_out) const
// {
//     sht_.forward( f_in.begin(), work_arr, 3*params_.n_surfs_, shc);
//     sht_.backward_du( shc, work_arr, 3*params_.n_surfs_, V1.begin());
//     sht_.backward_dv( shc, work_arr, 3*params_.n_surfs_, V2.begin());
    
//     GeometricDot(V1,cu_, S6);
//     GeometricDot(V2,cv_, div_f_out);
//     axpy((value_type)1.0,div_f_out,S6, div_f_out);

//     sht_.Filter(div_f_out.begin(),work_arr, params_.n_surfs_, shc, div_f_out.begin());
// }

// template <typename ScalarContainer, typename VectorContainer > 
// void Surface<ScalarContainer, VectorContainer>::Reparam()
// {
//     int iter = 0;
//     T vel = 2*params_.rep_max_vel_;

//     //up-sample
//     x_.getDevice().InterpSh(params_.p_, 3*params_.n_surfs_, x_.begin(), work_arr, shc, params_.rep_up_freq_, V10.begin());

//     while(iter++ < params_.rep_iter_max_ && vel > params_.rep_max_vel_ * params_.rep_max_vel_)
//     {
//         UpdateNormal();///@todo shc may have changed
//         x_.getDevice().Filter(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), alpha_q, work_arr, shc, V11.begin());

//         axpy((T) -1.0, V10, V11, V11);
//         GeometricDot(V13, V11, S10);
//         axpb((T) -1.0,S10, (T) 0.0, S10);
//         xvpw(S10, V13, V11, V11);//The correction velocity
//         axpy(params_.rep_ts_, V11, V10, V10);

//         GeometricDot(V11,V11,S10);
//         ///@todo The error should be relative
//         vel = S10.Max();
// #ifndef NDEBUG
//         cout<<" Reparam iteration "<<iter<<", max vel: "<<vel<<endl;
// #endif
//     }
    
//     ///inteprsh is buggy
//     x_.getDevice().InterpSh(params_.rep_up_freq_, 3*params_.n_surfs_, V10.begin(), work_arr, shc, params_.p_, V11.begin());
//     x_.getDevice().Memcpy(x_.begin(), V11.begin(), x_.getDataLength(), MemcpyDeviceToDevice);
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::UpdateNormal()
// {
//     x_.getDevice().ShAna(V10.begin(), work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, shc);
//     x_.getDevice().ShSynDu(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V11.begin());
//     x_.getDevice().ShSynDv(shc, work_arr, params_.rep_up_freq_, 3*params_.n_surfs_, V12.begin());

//     GeometricCross(V11,V12,V13);
//     GeometricDot(V13, V13, S10);
//     S10.Sqrt();
//     uyInv(V13,S10,V13);
// }

// template< typename ScalarContainer, typename VectorContainer >
// T Surface<ScalarContainer, VectorContainer>::Area()
// {
//     //Reduce(w_.begin(), quad_weights_, S1.getFunLength(), params_.n_surfs_, work_arr);
//     //return Max(work_arr, params_.n_surfs_);
//     return 0;
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::Volume()
// {
//     GeometricDot(x_,normal_,S1);
//     axpy((T) 1/3,S1, S1);
//     Reduce(S1.begin(), w_.begin(), quad_weights_, S1.getFunLength(), params_.n_surfs_, work_arr);
//     //     for(int ii=0;ii<params_.n_surfs_;++ii)
//     //         cout<<work_arr[ii]<<endl;
// }

// template< typename ScalarContainer, typename VectorContainer >
// void Surface<ScalarContainer, VectorContainer>::Populate(const T *centers)
// {
//     int length = this->x_.getFunLength();
//     for(int ii=1;ii<params_.n_surfs_;ii++)
//     {
//         int idx = 3 * ii * length;
//         x_.getDevice().Memcpy(x_.begin() + idx, x_.begin(), 3 * length, MemcpyDeviceToDevice);
//         x_.getDevice().axpb((T) 1.0, x_.begin() + idx                  , centers[3*ii  ], length, 1, x_.begin() + idx                  );
//         x_.getDevice().axpb((T) 1.0, x_.begin() + idx + length         , centers[3*ii+1], length, 1, x_.begin() + idx + length         );
//         x_.getDevice().axpb((T) 1.0, x_.begin() + idx + length + length, centers[3*ii+2], length, 1, x_.begin() + idx + length + length);
//     }

//     //treating the first surface
//     x_.getDevice().axpb((T) 1.0, x_.begin()                  , centers[0], length, 1, x_.begin()                  );
//     x_.getDevice().axpb((T) 1.0, x_.begin() + length         , centers[1], length, 1, x_.begin() + length         );
//     x_.getDevice().axpb((T) 1.0, x_.begin() + length + length, centers[2], length, 1, x_.begin() + length + length);
// }

// template<typename ScalarContainer, typename VectorContainer >
// T* Surface<ScalarContainer, VectorContainer>::GetCenters(T* cnts)
// {
//     UpdateFirstForms();
    
//     GeometricDot(x_, x_, S1);
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

//     GeometricDot(x_,normal_,S1);
//     axpb((T) 1/3,S1,(T) 0.0, S1);
//     x_.GetDevice().Reduce(S1.begin(), w_.begin(), quad_weights_, sc_len, nv, work_arr);
    
//     x_.GetDevice().Memcpy(work_arr + nv     , work_arr, nv, MemcpyDeviceToDevice); 
//     x_.GetDevice().Memcpy(work_arr + nv + nv, work_arr, nv, MemcpyDeviceTDevice); 
//     x_.GetDevice().ShufflePoints(work_arr, AxisMajor, nv, 1, work_arr + 3*nv);
//     x_.GetDevice().xyInv(cnts, work_arr + 3*nv, 3, nv, cnts);

//     return cnts;
// }
