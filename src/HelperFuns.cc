template<typename Container>
inline bool AreCompatible(const Container &lhs, const Container &rhs)
{
    return(lhs.GetFunLength() == rhs.GetFunLength()
        && lhs.GetNumFuns() == rhs.GetNumFuns());
}

template<typename ScalarContainer, typename VectorContainer>
inline bool AreCompatible(const VectorContainer &lhs, const  ScalarContainer &rhs)
{
    return(lhs.GetDevice() == rhs.GetDevice()
        && lhs.GetFunLength() == rhs.GetFunLength()
        && lhs.GetNumVecs() == rhs.GetNumFuns());
}

template<typename ScalarContainer, typename VectorContainer>
inline bool AreCompatible(const ScalarContainer &lhs, const VectorContainer &rhs)
{
    return(AreCompatible(rhs, lhs));
}

template<typename ScalarContainer>
void Sqrt(const ScalarContainer &x_in, ScalarContainer &sqrt_out)
{
    assert(AreCompatible(x_in, sqrt_out));

    x_in.GetDevice().Sqrt(x_in.begin(), x_in.Size(), sqrt_out.begin());
}

template<typename ScalarContainer>
void xy(const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &xy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xy_out));
    
    x_in.GetDevice().xy(x_in.begin(), y_in.begin(), x_in.Size(), xy_out.begin());
}

template<typename ScalarContainer>
void xyInv(const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &xyInv_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xyInv_out));
    
    x_in.GetDevice().xyInv(x_in.begin(), y_in.begin(), x_in.Size(), 
        xyInv_out.begin());
}

template<typename ScalarContainer>
void xInv(const ScalarContainer &x_in, ScalarContainer &xInv_out)
{
    assert(AreCompatible(x_in,xInv_out));

    x_in.GetDevice().xyInv(0, x_in.begin(), x_in.Size(), xInv_out.begin());
}

template<typename ScalarContainer>
void axpy(typename ScalarContainer::value_type a_in, 
    const ScalarContainer &x_in, const ScalarContainer &y_in, 
    ScalarContainer &axpy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,axpy_out));

    x_in.GetDevice().axpy(a_in, x_in.begin(), y_in.begin(), x_in.Size(), 
        axpy_out.begin());
}

template<typename ScalarContainer>
void axpy(typename ScalarContainer::value_type a_in, 
    const ScalarContainer &x_in, ScalarContainer &axpy_out)
{
    assert(AreCompatible(x_in,axpy_out));

    x_in.GetDevice().axpy<typename ScalarContainer::value_type>(
        a_in, x_in.begin(), 0, x_in.Size(), axpy_out.begin());
}

template<typename ScalarContainer>
void Reduce(const ScalarContainer &x_in, const ScalarContainer &w_in, 
    const ScalarContainer &quad_w_in, ScalarContainer &int_x_dw)
{
    assert(AreCompatible(x_in,w_in));
    assert(AreCompatible(w_in,int_x_dw));
    assert(quad_w_in.GetFunLength() == w_in.GetFunLength());
    assert(quad_w_in_in.GetNumFuns() >= 1);
    
    x_in.GetDevice().Reduce(x_in.begin(), w_in.begin(), quad_w_in.begin(), 
        x_in.GetFunLength(), x_in.GetNumFuns(), int_x_dw.begin());
}
 
template<typename ScalarContainer>
void Reduce(const ScalarContainer &w_in, const ScalarContainer &quad_w_in, 
    ScalarContainer &int_x_dw)
{
    assert(AreCompatible(w_in,int_x_dw));
    assert(quad_w_in.GetFunLength() == w_in.GetFunLength());
    assert(quad_w_in_in.GetNumFuns() >= 1);

    w_in.GetDevice().Reduce(0, w_in.begin(), quad_w_in.begin(), 
        w_in.GetFunLength(), w_in.GetNumFuns(), int_x_dw.begin());
}
   
template<typename ScalarContainer>
typename ScalarContainer::value_type Max(const ScalarContainer &x_in)
{
    return(x_in.GetDevice().Max(x_in.begin(), 
            x_in.GetFunLength() * x_in.GetNumFuns()));
}

template<typename ScalarContainer, typename VectorContainer>
inline void DotProduct(const VectorContainer &u_in,
    const VectorContainer &v_in, ScalarContainer &x_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,x_out));

    u_in.GetDevice().DotProduct(u_in.begin(), v_in.begin(),
        u_in.GetFunLength(), u_in.GetNumVecs(), x_out.begin());
}

template<typename VectorContainer>
inline void CrossProduct(const VectorContainer &u_in, 
    const VectorContainer &v_in, VectorContainer &w_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,w_out));

    u_in.GetDevice().CrossProduct(u_in.begin(), v_in.begin(), 
        u_in.GetFunLength(), u_in.GetNumVecs(), w_out.begin());
}

template<typename ScalarContainer, typename VectorContainer>
inline void uyInv(const VectorContainer &u_in, 
    const ScalarContainer &y_in, VectorContainer &uyInv_out)
{
    assert(AreCompatible(u_in,y_in));
    assert(AreCompatible(y_in,uyInv_out));

    u_in.GetDevice().uyInv(u_in.begin(), y_in.begin(), 
        u_in.GetFunLength(), u_in.GetNumVecs(), uyInv_out.begin());
}
    
template<typename ScalarContainer, typename VectorContainer>
inline void avpw(const ScalarContainer &a_in, 
    const VectorContainer &v_in,  const VectorContainer &w_in, 
    VectorContainer &avpw_out)
{
    assert(a_in.GetNumFuns() == v_in.GetNumVecs());
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,avpw_out));

    std::cout<<"THIS NEEDS TO BE IMPLEMENTED"<<std::endl;
    //u_in.GetDevice().avpw(a_in.begin(), v_in.begin(), 
    //    w_in.begin(), v_in.GetFunLength(), u_in.GetNumVecs(), avpw_out.begin());
}
 
template<typename ScalarContainer, typename VectorContainer>
inline void xvpw(const ScalarContainer &x_in, 
    const VectorContainer &v_in, const VectorContainer &w_in, 
    VectorContainer &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,xvpw_out));

    x_in.GetDevice().xvpw(x_in.begin(), v_in.begin(), 
        w_in.begin(), v_in.GetFunLength(), v_in.GetNumVecs(), xvpw_out.begin());
}

template<typename ScalarContainer, typename VectorContainer>
inline void xv(const ScalarContainer &x_in, 
    const VectorContainer &v_in, VectorContainer &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,xvpw_out));

    x_in.GetDevice().xvpw<typename ScalarContainer::value_type>(
        x_in.begin(), v_in.begin(), 0, v_in.GetFunLength(), 
        v_in.GetNumVecs(), xvpw_out.begin());
}
    
template<typename VectorContainer>
inline void ShufflePoints(const VectorContainer &u_in, 
    VectorContainer &u_out)
{
    assert(AreCompatible(u_in,u_out));
    size_t stride = u_in.GetFunLength();

    int dim = u_in.the_dim_;

    size_t dim1 = (u_in.GetOrder() == AxisMajor) ? dim : stride;
    size_t dim2 = (u_in.GetOrder() == AxisMajor) ? stride : dim;

    for(size_t ss=0; ss<u_in.GetNumVecs(); ++ss)
        u_in.GetDevice().transpose(u_in.begin() + dim * stride *ss
            , dim1, dim2, u_out.begin() + dim * stride * ss);
    
    u_out.SetOrder((u_in.GetOrder() == AxisMajor) ? PointMajor : AxisMajor);
}

///@todo this need to be implemented for the GPU
template<typename ScalarContainer>
typename ScalarContainer::value_type Dot(const ScalarContainer &x, 
    const ScalarContainer &y)
{
    typename ScalarContainer::value_type dot(0);
    for(size_t idx=0;idx<x.Size(); ++idx)
        dot+= x[idx] * y[idx];
    return dot;
}

///@todo This need a new data structure for the matrix, etc
template<typename ScalarContainer>
inline void CircShift(const typename ScalarContainer::value_type *x_in,
    int shift, ScalarContainer &x_out)
{
    int vec_length = x_out.GetFunLength();
    int n_vecs = x_out.GetNumFuns();
    shift = shift % vec_length;
    shift += (shift < 0) ?  vec_length : 0;
    
    int base_in, base_out;
    for (int ii = 0; ii < n_vecs; ii++) {
        base_out = ii * vec_length;
        base_in = base_out + vec_length - shift;
        x_out.GetDevice().Memcpy((float*) x_out.begin() + base_out, 
            x_in + base_in, sizeof(typename ScalarContainer::value_type) * shift, MemcpyDeviceToDevice);
        base_in = base_out;
        base_out += shift;
        x_out.GetDevice().Memcpy((float*) x_out.begin() + base_out, x_in + base_in, 
            sizeof(typename ScalarContainer::value_type) * (vec_length - shift), MemcpyDeviceToDevice);
    }
}
