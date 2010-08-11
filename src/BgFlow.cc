// Shear flow /////////////////////////////////////////////////////////////////
template<typename VecContainer>
ShearFlow<VecContainer>::ShearFlow(value_type shear_rate) :
    shear_rate_(shear_rate) {}


template<typename VecContainer>
void ShearFlow<VecContainer>::operator()(const VecContainer &pos, const value_type time,
        VecContainer &vel_inf)
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;
    
    axpy(0.0, pos, vel_inf);
    for(int ii=0;ii<n_surfs;ii++)
    {
        idx = pos.getTheDim() * ii * stride;
        pos.getDevice().Memcpy(vel_inf.begin() + idx, 
            pos.begin() + idx + stride + stride
            ,stride * sizeof(typename VecContainer::value_type), 
            MemcpyDeviceToDevice);
    }
    axpy(shear_rate_, vel_inf, vel_inf);
}

// Taylor vortex ////////////////////////////////////////////////////////////////
/**
 * @bug The data is manipulated directly that causes segmentation
 * fault on any device other than CPU.
 */
template<typename VecContainer>
TaylorVortex<VecContainer>::TaylorVortex(value_type strength, value_type x_period,
    value_type y_period) : 
    strength_(strength), x_period_(x_period), y_period_(y_period)
{
    assert( VecContainer::getDeviceType() == CPU );
}

template<typename VecContainer>
void TaylorVortex<VecContainer>::operator()(const VecContainer &pos, const value_type time,
    VecContainer &vel_inf)
{
    int n_surfs = pos.getNumSubs();
    int stride = pos.getStride();
    int idx;
    
    wrk_vec1_.replicate(pos);
    wrk_vec2_.replicate(pos);

    axpy( 2 * M_PI / x_period_, pos, wrk_vec1_);
    axpy( 2 * M_PI / y_period_, pos, wrk_vec2_);

    for ( size_t ii=0;ii<pos.size(); ++ii )
    {
        wrk_vec1_.begin()[ii] = cos(wrk_vec1_.begin()[ii]);
        wrk_vec2_.begin()[ii] = sin(wrk_vec2_.begin()[ii]); 
    }

    axpy(0.0, pos, vel_inf);
    for ( int ss=0; ss<n_surfs; ++ss )
        for ( int ii=0;ii<stride; ++ii)
        {
            idx = ss * DIM * stride + ii;
            vel_inf.begin()[idx          ] = - wrk_vec1_.begin()[idx] * wrk_vec2_.begin()[idx + stride];
            vel_inf.begin()[idx + stride ] = - wrk_vec2_.begin()[idx] * wrk_vec1_.begin()[idx + stride];

        }
    axpy(strength_, vel_inf, vel_inf); 
}
    

    

