/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */
template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::Surface(const Vec& x_in) :
    //upsample_freq_(2 * x_in.getShOrder()),
    rep_filter_freq_(x_in.getShOrder()/3),
    sht_(x_in.getShOrder()), ///@todo make sht_ autonomous
    //sht_upsample_(upsample_freq_),
    sht_rep_filter_(x_in.getShOrder(), rep_filter_freq_),
    containers_are_stale_(true),
    first_forms_are_stale_(true),
    second_forms_are_stale_(true),
    checked_out_work_sca_(0),
    checked_out_work_vec_(0)
{
    setPosition(x_in);
}

template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::~Surface()
{
    assert(!checked_out_work_sca_);
    assert(!checked_out_work_vec_);

    purgeTheWorkSpace();
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::setPosition(const Vec& x_in)
{
    containers_are_stale_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    x_.replicate(x_in);
    axpy(static_cast<value_type>(1), x_in, x_);
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer& Surface<ScalarContainer, 
                         VectorContainer>::getPositionModifiable()
{
    containers_are_stale_ = true;
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
void Surface<ScalarContainer, VectorContainer>::getSmoothedShapePosition(
    Vec &smthd_pos)
{
    Vec* wrk(produceVec(x_));
    Vec* shc(produceVec(x_));

    sht_rep_filter_.Filter(x_, *wrk, *shc, smthd_pos);

    recycleVec(wrk);
    recycleVec(shc);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::mapToTangentSpace(
    Vec &vec_fld)
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    Sca* scp(produceSca(x_));
    GeometricDot(vec_fld, normal_, *scp);
    axpy(static_cast<value_type>(-1.0), *scp, *scp);
    xvpw(*scp, normal_, vec_fld, vec_fld);
    
    recycleSca(scp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateFirstForms() const
{
    if(containers_are_stale_)
        checkContainers();
        
    Vec* wrk(produceVec(x_));
    Vec* shc(produceVec(x_));
    Vec* dif(produceVec(x_));
    Sca* scp(produceSca(x_));
    
    // Spherical harmonic coefficient
    sht_.FirstDerivatives(x_, *wrk, *shc, *dif, normal_);
    
    // First fundamental coefficients
    GeometricDot(*dif, *dif, E);
    GeometricDot(*dif, normal_, F);
    GeometricDot(normal_, normal_, G);
    
    // Area element
    xy(E, G, w_);
    xy(F, F, *scp);
    axpy(static_cast<value_type>(-1), *scp, w_, w_);

    // Dividing EFG by W^2
    xyInv(E, w_, E);
    xyInv(F, w_, F);
    xyInv(G, w_, G);   
    Sqrt(w_, w_);

    //Div and Grad coefficients
    xv(F, normal_, cu_);
    axpy(static_cast<value_type>(-1), cu_, cu_);
    xvpw(G, *dif, cu_, cu_);
    
    xv(F, *dif, cv_);
    axpy(static_cast<value_type>(-1), cv_, cv_);
    xvpw(E, normal_, cv_, cv_);

    GeometricCross(*dif, normal_, normal_);
    uyInv(normal_, w_, normal_);

    recycleVec(wrk);
    recycleVec(shc);
    recycleVec(dif);
    recycleSca(scp);    

    first_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::updateAll() const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Vec* wrk(produceVec(x_));
    Vec* shc(produceVec(x_));
    Vec* dif(produceVec(x_));

    sht_.forward(x_, *wrk, *shc);
    
    sht_.backward_duv(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, h_);

    xy(h_, h_, k_);
    axpy(static_cast<value_type>(-1), k_, k_);

    xy(F, h_, h_);
    axpy(static_cast<value_type>(-1), h_, h_);

    Sca* L(produceSca(x_));
    sht_.backward_d2u(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *L);
    
    Sca* N(produceSca(x_));
    xy(G, *L, *N);
    axpy(static_cast<value_type>(.5), *N, h_, h_); 

    
    sht_.backward_d2v(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *N);

    xy(*L, *N, *L);
    axpy(static_cast<value_type>(1), *L, k_, k_); 

    xyInv(k_,w_,k_);
    xyInv(k_,w_,k_);

    xy(E, *N, *N);     
    axpy(static_cast<value_type>(.5), *N, h_, h_); 

    recycleSca(L);
    recycleSca(N);

    sht_.Filter(k_, *wrk, *shc, k_);
    sht_.Filter(h_, *wrk, *shc, h_);

    recycleVec(wrk);
    recycleVec(shc);
    recycleVec(dif);

    second_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::grad(const ScalarContainer 
    &f_in, VectorContainer &grad_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca* scw1(produceSca(x_));
    Sca* scw2(produceSca(x_));
    Vec* shc(produceVec(x_));
    Vec* wrk(produceVec(x_));

    sht_.FirstDerivatives(f_in, *wrk, *shc, *scw1, *scw2);

    xv(*scw1, cu_, grad_f_out);
    recycleSca(scw1);

    xvpw(*scw2, cv_, grad_f_out, grad_f_out);
    recycleSca(scw2);

    sht_.Filter(grad_f_out, *wrk, *shc, grad_f_out);
    recycleVec(shc);
    recycleVec(wrk);
}

template <typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::div(const VectorContainer 
    &f_in, ScalarContainer &div_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Vec* dif(produceVec(x_));
    Vec* shc(produceVec(x_));
    Vec* wrk(produceVec(x_));
    Sca* scw(produceSca(x_));
    
    sht_.forward(f_in, *wrk, *shc);
    sht_.backward_du(*shc, *wrk, *dif);
    GeometricDot(*dif, cu_, *scw);
    
    sht_.backward_dv(*shc, *wrk, *dif);
    GeometricDot(*dif, cv_, div_f_out);
    recycleVec(dif);

    axpy(static_cast<value_type>(1),div_f_out, *scw, div_f_out);
    recycleSca(scw);

    sht_.Filter(div_f_out, *wrk, *shc, div_f_out);
    recycleVec(shc);
    recycleVec(wrk);  
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::area(ScalarContainer 
    &area_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    integrator_(w_, area_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::volume(ScalarContainer 
    &vol_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    Sca* scw(produceSca(x_));
    GeometricDot(x_,normal_,*scw);
    axpy(static_cast<value_type>(1)/3, *scw, *scw);

    integrator_(*scw, w_, vol_out);
    recycleSca(scw);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::getCenters(Vec &centers) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca* scw(produceSca(x_));
    GeometricDot(x_, x_, *scw);
    axpy(static_cast<value_type>(.5), *scw, *scw);
    
    Vec* vcw(produceVec(x_));
    xv(*scw, normal_, *vcw);
    recycleSca(scw);

    integrator_(*vcw, w_, centers);
    recycleVec(vcw);

    scw = produceSca(centers);
    volume(*scw);
    uyInv(centers, *scw, centers);
    recycleSca(scw);
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
    
    E.replicate(x_);    
    F.replicate(x_);    
    G.replicate(x_);    
    
    containers_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
ScalarContainer* Surface<ScalarContainer, VectorContainer>::produceSca(
    const VectorContainer &ref) const
{
    Sca* scp;
    
    if(scalar_work_q_.empty())
        scp = new ScalarContainer;
    else
    {
        scp = scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    scp->replicate(ref);
    ++checked_out_work_sca_;
    return(scp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::recycleSca(
    ScalarContainer* scp) const
{
    scalar_work_q_.push(scp);
    --checked_out_work_sca_;
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer* Surface<ScalarContainer, VectorContainer>::produceVec(
    const VectorContainer &ref) const
{
    Vec* vcp;
    
    if(vector_work_q_.empty())
        vcp = new VectorContainer;
    else
    {
        vcp = vector_work_q_.front();
        vector_work_q_.pop();
    }
    
    vcp->replicate(ref);
    ++checked_out_work_vec_;
    
    return(vcp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::recycleVec(
    VectorContainer* vcp) const
{
    vector_work_q_.push(vcp);
    --checked_out_work_vec_;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::purgeTheWorkSpace() const
{
    while ( !scalar_work_q_.empty() )
    {
        delete scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    while ( !vector_work_q_.empty() )
    {
        delete vector_work_q_.front();
        vector_work_q_.pop();
    }
}
