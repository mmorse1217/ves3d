/**
 * @file   Surface.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:37:21 2010
 * 
 * @brief  The implementation of the surface class.
 */
template <typename ScalarContainer, typename VectorContainer>  
Surface<ScalarContainer, VectorContainer>::Surface(
    const Vec_t& x_in, OperatorsMats<value_type, device_type> &mats) :
    upsample_freq_(2 * x_in.getShOrder()),
    rep_filter_freq_(x_in.getShOrder()/3),
    sht_(x_in.getShOrder(), mats.mats_p_), ///@todo make sht_ autonomous
    sht_rep_filter_(x_in.getShOrder(), mats.mats_p_, rep_filter_freq_),
    sht_rep_upsample_(upsample_freq_, mats.mats_p_up_),
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
void Surface<ScalarContainer, VectorContainer>::setPosition(const Vec_t& x_in)
{
    containers_are_stale_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    x_.replicate(x_in);
    axpy(static_cast<value_type>(1), x_in, x_);
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer& Surface<ScalarContainer, VectorContainer>::getPositionModifiable()
{
    containers_are_stale_ = true;
    first_forms_are_stale_ = true;
    second_forms_are_stale_ = true;

    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, VectorContainer>::getPosition() const
{
    return(x_);
}

template <typename ScalarContainer, typename VectorContainer>  
const VectorContainer& Surface<ScalarContainer, VectorContainer>::getNormal() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(normal_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer,VectorContainer>::getAreaElement() const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    return(w_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer, VectorContainer>::getMeanCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(h_);
}

template <typename ScalarContainer, typename VectorContainer>  
const ScalarContainer& Surface<ScalarContainer,VectorContainer>::getGaussianCurv() const
{
    if(second_forms_are_stale_)
        updateAll();
    return(k_);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
getSmoothedShapePosition(Vec_t &smthd_pos) const
{
    Vec_t* wrk(checkoutVec());
    Vec_t* shc(checkoutVec());
    
    sht_rep_filter_.lowPassFilter(x_, *wrk, *shc, smthd_pos);
    
    recycleVec(wrk);
    recycleVec(shc);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
mapToTangentSpace(Vec_t &vec_fld) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca_t* scp(checkoutSca());
    Vec_t* wrk(checkoutVec());
    Vec_t* shc(checkoutVec());
    Vec_t* fld(checkoutVec());

    //up-sampling
    int usf(sht_rep_upsample_.getShOrder());

    scp->resize(scp->getNumSubs(), usf);
    wrk->resize(wrk->getNumSubs(), usf);
    shc->resize(shc->getNumSubs(), usf);
    fld->resize(fld->getNumSubs(), usf);
    normal_.resize(normal_.getNumSubs(), usf);

    Resample(vec_fld, sht_, sht_rep_upsample_, *shc, *wrk, *fld);
     Resample(normal_, sht_, sht_rep_upsample_, *shc, *wrk, normal_);
    
    //re-normalizing
    GeometricDot(normal_, normal_, *scp);
    Sqrt(*scp, *scp);
    uyInv(normal_, *scp, normal_);
    GeometricDot(normal_, normal_, *scp);
    
    //mapping to tangent
    GeometricDot(*fld, normal_, *scp);
    axpy(static_cast<value_type>(-1.0), *scp, *scp);
    xvpw(*scp, normal_, *fld, *fld);
    
    //down-sampling
    Resample(*fld   , sht_rep_upsample_, sht_, *shc, *wrk, vec_fld);
    Resample(normal_, sht_rep_upsample_, sht_, *shc, *wrk, normal_);
    normal_.resize(normal_.getNumSubs(), sht_.getShOrder());

    recycleSca(scp);
    recycleVec(wrk);
    recycleVec(shc);
    recycleVec(fld);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
updateFirstForms() const
{
    if(containers_are_stale_)
        checkContainers();
        
    Vec_t* wrk(checkoutVec());
    Vec_t* shc(checkoutVec());
    Vec_t* dif(checkoutVec());
    Sca_t* scp(checkoutSca());
    
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

    recycle(wrk);
    recycle(shc);
    recycle(dif);
    recycle(scp);    

    first_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
updateAll() const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Vec_t* wrk(checkoutVec());
    Vec_t* shc(checkoutVec());
    Vec_t* dif(checkoutVec());

    sht_.forward(x_, *wrk, *shc);
    
    sht_.backward_duv(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, h_);

    xy(h_, h_, k_);
    axpy(static_cast<value_type>(-1), k_, k_);

    xy(F, h_, h_);
    axpy(static_cast<value_type>(-1), h_, h_);

    Sca_t* L(checkoutSca());
    sht_.backward_d2u(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *L);
    
    Sca_t* N(checkoutSca());
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

    recycle(L);
    recycle(N);

    sht_.lowPassFilter(k_, *wrk, *shc, k_);
    sht_.lowPassFilter(h_, *wrk, *shc, h_);

    recycle(wrk);
    recycle(shc);
    recycle(dif);

    second_forms_are_stale_ = false;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
grad(const ScalarContainer &f_in, VectorContainer &grad_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca_t* scw1(checkoutSca());
    Sca_t* scw2(checkoutSca());
    Vec_t* shc(checkoutVec());
    Vec_t* wrk(checkoutVec());

    sht_.FirstDerivatives(f_in, *wrk, *shc, *scw1, *scw2);

    xv(*scw1, cu_, grad_f_out);
    recycle(scw1);

    xvpw(*scw2, cv_, grad_f_out, grad_f_out);
    recycle(scw2);

    sht_.lowPassFilter(grad_f_out, *wrk, *shc, grad_f_out);
    recycle(shc);
    recycle(wrk);
}

template <typename ScalarContainer, typename VectorContainer> 
void Surface<ScalarContainer, VectorContainer>::
div(const VectorContainer &f_in, ScalarContainer &div_f_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Vec_t* dif(checkoutVec());
    Vec_t* shc(checkoutVec());
    Vec_t* wrk(checkoutVec());
    Sca_t* scw(checkoutSca());
    
    sht_.forward(f_in, *wrk, *shc);
    sht_.backward_du(*shc, *wrk, *dif);
    GeometricDot(*dif, cu_, *scw);
    
    sht_.backward_dv(*shc, *wrk, *dif);
    GeometricDot(*dif, cv_, div_f_out);
    recycle(dif);

    axpy(static_cast<value_type>(1),div_f_out, *scw, div_f_out);
    recycle(scw);

    sht_.lowPassFilter(div_f_out, *wrk, *shc, div_f_out);
    recycle(shc);
    recycle(wrk);  
}

template< typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
area(ScalarContainer &area_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    integrator_(w_, area_out);
}

template< typename ScalarContainer, typename VectorContainer >
void Surface<ScalarContainer, VectorContainer>::
volume(ScalarContainer &vol_out) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    Sca_t* scw(checkoutSca());
    GeometricDot(x_,normal_,*scw);
    axpy(static_cast<value_type>(1)/3, *scw, *scw);

    integrator_(*scw, w_, vol_out);
    recycleSca(scw);
}

template< typename ScalarContainer, typename VectorContainer>
void Surface<ScalarContainer, VectorContainer>::
getCenters(Vec_t &centers) const
{
    if(first_forms_are_stale_)
        updateFirstForms();

    Sca_t* scw(checkoutSca());
    GeometricDot(x_, x_, *scw);
    axpy(static_cast<value_type>(.5), *scw, *scw);
    
    Vec_t* vcw(checkoutVec());
    xv(*scw, normal_, *vcw);
    recycleSca(scw);

    integrator_(*vcw, w_, centers);
    recycleVec(vcw);

    scw = checkoutSca();
    scw->replicate(centers);
    volume(*scw);
    uyInv(centers, *scw, centers);
    recycleSca(scw);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
checkContainers() const
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
ScalarContainer* Surface<ScalarContainer, VectorContainer>::
checkoutSca() const
{
    Sca_t* scp;
    
    if(scalar_work_q_.empty())
        scp = new ScalarContainer;
    else
    {
        scp = scalar_work_q_.front();
        scalar_work_q_.pop();
    }
    
    scp->replicate(x_);
    ++checked_out_work_sca_;
    return(scp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
recycle(Sca_t* scp) const
{
    scalar_work_q_.push(scp);
    --checked_out_work_sca_;
}

template <typename ScalarContainer, typename VectorContainer>  
VectorContainer* Surface<ScalarContainer, VectorContainer>::
checkoutVec() const
{
    Vec_t* vcp;
    
    if(vector_work_q_.empty())
        vcp = new VectorContainer;
    else
    {
        vcp = vector_work_q_.front();
        vector_work_q_.pop();
    }
    
    vcp->replicate(x_);
    ++checked_out_work_vec_;
    
    return(vcp);
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::recycle(
    VectorContainer* vcp) const
{
    vector_work_q_.push(vcp);
    --checked_out_work_vec_;
}

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
purgeTheWorkSpace() const
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

template <typename ScalarContainer, typename VectorContainer>  
void Surface<ScalarContainer, VectorContainer>::
linearizedMeanCurv(const Vec_t &x_new, Sca_t &h_lin) const
{
    if(first_forms_are_stale_)
        updateFirstForms();
    
    Vec_t* wrk(checkoutVec());
    Vec_t* shc(checkoutVec());
    Vec_t* dif(checkoutVec());
    Sca_t* scw(checkoutSca());

    sht_.forward(x_new, *wrk, *shc);   

    //duu
    sht_.backward_d2u(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, h_lin);
    xy(G, h_lin, h_lin);

    //duv
    sht_.backward_duv(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *scw);
    xy(F, *scw, *scw);
    axpy(static_cast<value_type>(-2), *scw, h_lin, h_lin);

    //dvv
    sht_.backward_d2v(*shc, *wrk, *dif);
    GeometricDot(*dif, normal_, *scw);
    xy(E, *scw, *scw);
    axpy(static_cast<value_type>(1), *scw, h_lin, h_lin);

    axpy(static_cast<value_type>(.5), h_lin, h_lin);

    sht_.lowPassFilter(h_lin, *wrk, *shc, h_lin);

    recycle(wrk);
    recycle(shc);
    recycle(dif);
    recycle(scw);
}
