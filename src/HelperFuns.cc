template<typename lhsContainer, typename rhsContainer>
inline bool AreCompatible(const lhsContainer &lhs, const  rhsContainer &rhs)
{
    return(lhs.getDevice() == rhs.getDevice()
        && lhs.getStride() == rhs.getStride()
        && lhs.getNumSubs() == rhs.getNumSubs());
}

template<typename ScalarContainer>
void Sqrt(const ScalarContainer &x_in, ScalarContainer &sqrt_out)
{
    assert(AreCompatible(x_in, sqrt_out));

    x_in.getDevice().Sqrt(x_in.begin(), x_in.size(), sqrt_out.begin());
}

template<typename ScalarContainer>
void xy(const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &xy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xy_out));
    
    x_in.getDevice().xy(x_in.begin(), y_in.begin(), x_in.size(), xy_out.begin());
}

template<typename ScalarContainer>
void xyInv(const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &xyInv_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xyInv_out));
    
    x_in.getDevice().xyInv(x_in.begin(), y_in.begin(), x_in.size(), 
        xyInv_out.begin());
}

template<typename ScalarContainer>
void xInv(const ScalarContainer &x_in, ScalarContainer &xInv_out)
{
    assert(AreCompatible(x_in,xInv_out));

    x_in.getDevice().xyInv(0, x_in.begin(), x_in.size(), xInv_out.begin());
}

template<typename ScalarContainer>
void axpy(typename ScalarContainer::value_type a_in, 
    const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &axpy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,axpy_out));

    x_in.getDevice().axpy(a_in, x_in.begin(), y_in.begin(), x_in.size(), 
        axpy_out.begin());
}

template<typename ScalarContainer>
void axpy(typename ScalarContainer::value_type a_in, 
    const ScalarContainer &x_in, ScalarContainer &axpy_out)
{
    assert(AreCompatible(x_in,axpy_out));

    x_in.getDevice().axpy<typename ScalarContainer::value_type>(
        a_in, x_in.begin(), 0, x_in.size(), axpy_out.begin());
}


template<typename ScalarContainer, typename IntegrandContainer>
void Reduce(const IntegrandContainer &x_in, 
    const ScalarContainer &w_in, const ScalarContainer &quad_w_in, 
    IntegrandContainer &x_dw)
{
    assert(AreCompatible(x_in,w_in));
    assert(quad_w_in.getStride() == w_in.getStride());
    assert(quad_w_in.getNumSubs() >= 1);
    
    x_in.getDevice().Reduce(x_in.begin(), x_in.getTheDim(), w_in.begin(), 
        quad_w_in.begin(), x_in.getStride(), x_in.getNumSubs(), x_dw.begin());
}

template<typename Container>
void Reduce(const Container &w_in, const Container &quad_w_in, 
    Container &dw)
{
    assert(quad_w_in.getStride() == w_in.getStride());
    assert(quad_w_in.getNumSubs() >= 1);
    
    w_in.getDevice().Reduce(static_cast<typename Container::value_type* >(0), 
        0, w_in.begin(), quad_w_in.begin(), w_in.getStride(), 
        w_in.getNumSubs(), dw.begin());
}
 
template<typename ScalarContainer>
typename ScalarContainer::value_type Max(const ScalarContainer &x_in)
{
    return(x_in.getDevice().Max(x_in.begin(), 
            x_in.getStride() * x_in.getNumSubs()));
}

template<typename ScalarContainer, typename VectorContainer>
inline void GeometricDot(const VectorContainer &u_in,
    const VectorContainer &v_in, ScalarContainer &x_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,x_out));

    u_in.getDevice().DotProduct(u_in.begin(), v_in.begin(),
        u_in.getStride(), u_in.getNumSubs(), x_out.begin());
}

template<typename VectorContainer>
inline void GeometricCross(const VectorContainer &u_in, 
    const VectorContainer &v_in, VectorContainer &w_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,w_out));

    u_in.getDevice().CrossProduct(u_in.begin(), v_in.begin(), 
        u_in.getStride(), u_in.getNumSubs(), w_out.begin());
}

template<typename ScalarContainer, typename VectorContainer>
inline void uyInv(const VectorContainer &u_in, 
    const ScalarContainer &y_in, VectorContainer &uyInv_out)
{
    assert(AreCompatible(u_in, y_in));
    assert(AreCompatible(y_in,uyInv_out));

    u_in.getDevice().uyInv(u_in.begin(), y_in.begin(), 
        u_in.getStride(), u_in.getNumSubs(), uyInv_out.begin());
}
    
template<typename ScalarContainer, typename VectorContainer>
inline void avpw(const ScalarContainer &a_in, 
    const VectorContainer &v_in,  const VectorContainer &w_in, 
    VectorContainer &avpw_out)
{
    assert(a_in.getNumFuns() == v_in.getNumSubs());
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,avpw_out));

    std::cerr<<"THIS NEEDS TO BE IMPLEMENTED"<<std::endl;
    //u_in.getDevice().avpw(a_in.begin(), v_in.begin(), 
    //    w_in.begin(), v_in.getStride(), u_in.getNumSubs(), avpw_out.begin());
}
 
template<typename ScalarContainer, typename VectorContainer>
inline void xvpw(const ScalarContainer &x_in, 
    const VectorContainer &v_in, const VectorContainer &w_in, 
    VectorContainer &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,xvpw_out));

    x_in.getDevice().xvpw(x_in.begin(), v_in.begin(), 
        w_in.begin(), v_in.getStride(), v_in.getNumSubs(), xvpw_out.begin());
}

template<typename ScalarContainer, typename VectorContainer>
inline void xv(const ScalarContainer &x_in, 
    const VectorContainer &v_in, VectorContainer &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,xvpw_out));

    x_in.getDevice().xvpw<typename ScalarContainer::value_type>(
        x_in.begin(), v_in.begin(), 0, v_in.getStride(), 
        v_in.getNumSubs(), xvpw_out.begin());
}
    
template<typename VectorContainer>
inline void ShufflePoints(const VectorContainer &u_in, 
    VectorContainer &u_out)
{
    assert(AreCompatible(u_in,u_out));
    size_t stride = u_in.getStride();

    int dim = u_in.the_dim_;

    size_t dim1 = (u_in.getOrder() == AxisMajor) ? dim : stride;
    size_t dim2 = (u_in.getOrder() == AxisMajor) ? stride : dim;

    for(size_t ss=0; ss<u_in.getNumSubs(); ++ss)
        u_in.getDevice().transpose(u_in.begin() + dim * stride *ss
            , dim1, dim2, u_out.begin() + dim * stride * ss);
    
    u_out.SetOrder((u_in.getOrder() == AxisMajor) ? PointMajor : AxisMajor);
}

///@todo This need a new data structure for the matrix, etc
template<typename ScalarContainer>
inline void CircShift(const typename ScalarContainer::value_type *x_in,
    int shift, ScalarContainer &x_out)
{
    int vec_length = x_out.getStride();
    int n_vecs = x_out.getNumSubs();
    shift = shift % vec_length;
    shift += (shift < 0) ?  vec_length : 0;

    int base_in, base_out;
    for (int ii = 0; ii < n_vecs; ii++) {
        base_out = ii * vec_length;
        base_in = base_out + vec_length - shift;
        x_out.getDevice().Memcpy((float*) x_out.begin() + base_out, 
            x_in + base_in, sizeof(typename ScalarContainer::value_type) * shift, MemcpyDeviceToDevice);
        base_in = base_out;
        base_out += shift;
        x_out.getDevice().Memcpy((float*) x_out.begin() + base_out, x_in + base_in, 
            sizeof(typename ScalarContainer::value_type) * (vec_length - shift), MemcpyDeviceToDevice);
    }
}

///@todo this need to be implemented for the GPU
template<typename ScalarContainer>
typename ScalarContainer::value_type AlgebraicDot(const ScalarContainer &x, 
    const ScalarContainer &y)
{
    typename ScalarContainer::value_type dot(0);
    for(size_t idx=0;idx<x.size(); ++idx)
        dot+= x[idx] * y[idx];
    return dot;
}
