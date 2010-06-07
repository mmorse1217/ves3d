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
    sht_.backward_du(shc.begin(), work_arr.begin(), nv3, Xu.begin());
    sht_.backward_dv(shc.begin(), work_arr.begin(), nv3, Xv.begin());
    
    // First fundamental coefficients
    GeometricDot(Xu, Xu, E);
    GeometricDot(Xu, Xv, F);
    GeometricDot(Xv, Xv, G);
    GeometricCross(Xu, Xv, normal_);
    
    // Area element
    GeometricDot(normal_, normal_, w_);//w_ = W^2 for now
    
    // Dividing EFG by W^2
    xyInv(E, w_, E);
    xyInv(F, w_, F);
    xyInv(G, w_, G);
    Sqrt(w_, w_);
    uyInv(normal_, w_, normal_);

    //Div and Grad coefficients
    xv(F, Xv, cu_);
    axpy((value_type) -1.0,cu_, cu_);
    xvpw(G, Xu, cu_, cu_);
    
    xv(F, Xu, cv_);
    axpy((value_type) -1.0,cv_, cv_);
    xvpw(E, Xv, cv_, cv_);

    first_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateAll() const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    size_t nv3 = x_.getTheDim() * x_.getNumSubs();
    
    // Second derivatives
    sht_.forward(x_.begin(), work_arr.begin(), nv3, shc.begin());

    sht_.backward_d2u(shc.begin(), work_arr.begin(), nv3, Xu.begin());
    GeometricDot(Xu, normal_, L);
    
    sht_.backward_duv(shc.begin(), work_arr.begin(), nv3, Xu.begin());
    GeometricDot(Xu, normal_, M);

    sht_.backward_d2v(shc.begin(), work_arr.begin(), nv3, Xu.begin());
    GeometricDot(Xu, normal_, N);
    
    // Gaussian curvature
    xy(L,N,k_);
    xy(M,M,h_);
    axpy((value_type)-1.0,h_,k_,k_);
    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);
    
    sht_.Filter(k_.begin(),work_arr.begin(), k_.getNumSubs(), 
        shc.begin(), k_.begin());
    
    // Mean curvature
    xy(E,N,h_);
    axpy((value_type) .5, h_, h_);
    
    xy(F,M,N);
    axpy((value_type)-1.0, N, h_, h_);
    
    xy(G,L,N);
    axpy((value_type).5 , N, h_, h_);

    sht_.Filter(h_.begin(),work_arr.begin(), h_.getNumSubs(), 
        shc.begin(), h_.begin());

    second_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::grad(const ScalarContainer 
    &f_in, VectorContainer &grad_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    sht_.forward(f_in.begin(), work_arr.begin(), f_in.getNumSubs(), shc.begin());
    sht_.backward_du(shc.begin(), work_arr.begin(), f_in.getNumSubs(), M.begin());
    sht_.backward_dv(shc.begin(), work_arr.begin(), f_in.getNumSubs(), N.begin());

    xv(M, cu_, grad_f_out);
    xvpw(N, cv_, grad_f_out, grad_f_out);
    sht_.Filter(grad_f_out.begin(), work_arr.begin(), 
        grad_f_out.getTheDim() * grad_f_out.getNumSubs(),
        shc.begin(), grad_f_out.begin());
}

template <typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::div(const VectorContainer 
    &f_in, ScalarContainer &div_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    size_t nv3 = f_in.getTheDim() * f_in.getNumSubs();

    sht_.forward(f_in.begin(), work_arr.begin(), nv3, shc.begin());
    sht_.backward_du( shc.begin(), work_arr.begin(), nv3, Xu.begin());
    sht_.backward_dv( shc.begin(), work_arr.begin(), nv3, Xv.begin());
    
    GeometricDot(Xu,cu_, N);
    GeometricDot(Xv,cv_, div_f_out);
    axpy((value_type)1.0,div_f_out, N, div_f_out);

    sht_.Filter(div_f_out.begin(),work_arr.begin(), div_f_out.getNumSubs(), 
        shc.begin(), div_f_out.begin());
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::area(ScalarContainer &area_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    integrator_(w_, area_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::volume(ScalarContainer &vol_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    GeometricDot(x_,normal_,N);
    axpy(static_cast<typename ScalarContainer::value_type>(1)/3, N, N);

    integrator_(N, w_, vol_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::populate(const Sca &centers)
{
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::getCenters(Sca &centers) const
{
}
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

//     x_.GetDevice().Memcpy(work_arr.begin()            , w_.begin(), len, MemcpyDeviceToDevice);
//     x_.GetDevice().Memcpy(work_arr.begin() + len      , w_.begin(), len, MemcpyDeviceToDevice);
//     x_.GetDevice().Memcpy(work_arr.begin() + len + len, w_.begin(), len, MemcpyDeviceToDevice);

//     x_.GetDevice().ShufflePoints(work_arr.begin()                   , AxisMajor , len   , 1 , work_arr.begin() + len + len + len);
//     x_.GetDevice().ShufflePoints(work_arr.begin() + len + len + len , PointMajor, sc_len, nv, work_arr.begin());

//     x_.GetDevice().Reduce(NULL, work_arr.begin(), quad_weights_, sc_len, 3*nv, cnts);
//     x_.GetDevice().Reduce(V1.begin(), work_arr.begin(), quad_weights_, sc_len, 3*nv, cnts);

//     GeometricDot(x_,normal_,S1);
//     axpb((T) 1/3,S1,(T) 0.0, S1);
//     x_.GetDevice().Reduce(S1.begin(), w_.begin(), quad_weights_, sc_len, nv, work_arr.begin());
    
//     x_.GetDevice().Memcpy(work_arr.begin() + nv     , work_arr.begin(), nv, MemcpyDeviceToDevice); 
//     x_.GetDevice().Memcpy(work_arr.begin() + nv + nv, work_arr.begin(), nv, MemcpyDeviceTDevice); 
//     x_.GetDevice().ShufflePoints(work_arr.begin(), AxisMajor, nv, 1, work_arr.begin() + 3*nv);
//     x_.GetDevice().xyInv(cnts, work_arr.begin() + 3*nv, 3, nv, cnts);

//     return cnts;
// }

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
    E.replicate(x_);    
    F.replicate(x_);    
    G.replicate(x_);    
    L.replicate(x_);    
    M.replicate(x_);    
    N.replicate(x_);    
    Xu.replicate(x_);    
    Xv.replicate(x_);

    position_has_changed_outside_ = false;
}
