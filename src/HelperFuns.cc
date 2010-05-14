template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Scalars<T,DT> &lhs, const Scalars<T,DT> &rhs)
{
    return(lhs.GetDevicePtr()  == rhs.GetDevicePtr()
        && lhs.GetFunLength()     == rhs.GetFunLength()
        && lhs.GetNumFuns()    == rhs.GetNumFuns());
}

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Vectors<T,DT> &lhs, const  Scalars<T,DT> &rhs)
{
    return(lhs.GetDevicePtr()  == rhs.GetDevicePtr()
        && lhs.GetFunLength()     == rhs.GetFunLength()
        && lhs.GetNumVecs()    == rhs.GetNumFuns());
}

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Scalars<T,DT> &lhs, const Vectors<T,DT> &rhs)
{
    return(AreCompatible(rhs, lhs));
}

template<typename T, enum DeviceType DT>
inline bool AreCompatible(const Vectors<T,DT> &lhs, const Vectors<T,DT> &rhs)
{
    return(lhs.GetDevicePtr()  == rhs.GetDevicePtr()
        && lhs.GetFunLength()     == rhs.GetFunLength()
        && lhs.GetNumFuns()    == rhs.GetNumFuns());
}

template<typename T, enum DeviceType DT>
void Sqrt(const Scalars<T,DT> &x_in, Scalars<T,DT> &sqrt_out)
{
    assert(AreCompatible(x_in,sqrt_out));

    x_in.GetDevicePtr()->Sqrt(x_in.begin(), x_in.Size(),
        sqrt_out.begin());
}

template<typename T, enum DeviceType DT>
void xy(const Scalars<T,DT> &x_in, 
        const Scalars<T,DT> &y_in, 
              Scalars<T,DT> &xy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xy_out));

    x_in.GetDevicePtr()->xy(x_in.begin(), y_in.begin(), 
        x_in.Size(), xy_out.begin());
}

template<typename T, enum DeviceType DT>
void xyInv(const Scalars<T,DT> &x_in, 
           const Scalars<T,DT> &y_in, 
                 Scalars<T,DT> &xyInv_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,xyInv_out));

    x_in.GetDevicePtr()->xyInv(x_in.begin(), y_in.begin(), 
        x_in.Size(), xyInv_out.begin());
}

template<typename T, enum DeviceType DT>
void xInv(const Scalars<T,DT> &x_in, 
                Scalars<T,DT> &xInv_out)
{
    assert(AreCompatible(x_in,xInv_out));

    x_in.GetDevicePtr()->xyInv(NULL, x_in.begin(), 
        x_in.Size(), xInv_out.begin());
}

template<typename T, enum DeviceType DT>
void axpy(T a_in, const Scalars<T,DT> &x_in, 
                  const Scalars<T,DT> &y_in, 
                        Scalars<T,DT> &axpy_out)
{
    assert(AreCompatible(x_in,y_in));
    assert(AreCompatible(y_in,axpy_out));

    x_in.GetDevicePtr()->axpy(a_in, x_in.begin(), y_in.begin(), 
        x_in.Size(), axpy_out.begin());
}

template<typename T, enum DeviceType DT>
void axpy(T a_in, const Scalars<T,DT> &x_in, Scalars<T,DT> &axpy_out)
{
    assert(AreCompatible(x_in,axpy_out));

    x_in.GetDevicePtr()->axpy<T>(a_in, x_in.begin(), NULL, 
        x_in.Size(), axpy_out.begin());
}

template<typename T, enum DeviceType DT>
void Reduce(const Scalars<T,DT> &x_in, 
            const Scalars<T,DT> &w_in, 
            const Scalars<T,DT> &quad_w_in, 
                  Scalars<T,DT> &int_x_dw)
{
    assert(AreCompatible(x_in,w_in));
    assert( AreCompatible(w_in,int_x_dw));
    assert(quad_w_in.GetFunLength() == w_in.GetFunLength());
    assert(quad_w_in_in.GetNumFuns() >= 1);
    
    x_in.GetDevicePtr()->Reduce(x_in.begin(), w_in.begin(), quad_w_in.begin(), 
        x_in.GetFunLength(), x_in.GetNumFuns(), int_x_dw.begin());
}
 
template<typename T, enum DeviceType DT>
void Reduce(const Scalars<T,DT> &w_in, 
            const Scalars<T,DT> &quad_w_in, 
                  Scalars<T,DT> &int_x_dw)
{
    assert(AreCompatible(w_in,int_x_dw));
    assert(quad_w_in.GetFunLength() == w_in.GetFunLength());
    assert(quad_w_in_in.GetNumFuns() >= 1);

    w_in.GetDevicePtr()->Reduce(NULL, w_in.begin(), quad_w_in.begin(), 
        w_in.GetFunLength(), w_in.GetNumFuns(), int_x_dw.begin());
}
   
template<typename T, enum DeviceType DT>
T Max(const Scalars<T,DT> &x_in)
{
    return(x_in.GetDevicePtr()->Max(x_in.begin(), 
            x_in.GetFunLength() * x_in.GetNumFuns()));
}
 
template<typename T, enum DeviceType DT>
inline void DotProduct(const Vectors<T,DT> &u_in, 
                       const Vectors<T,DT> &v_in, 
                             Scalars<T,DT> &x_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,x_out));

    u_in.GetDevicePtr()->DotProduct(u_in.begin(), v_in.begin(),
        u_in.GetFunLength(), u_in.GetNumVecs(), x_out.begin());
}

template<typename T, enum DeviceType DT>
inline void CrossProduct(const Vectors<T,DT> &u_in, 
                         const Vectors<T,DT> &v_in, 
                               Vectors<T,DT> &w_out)
{
    assert(AreCompatible(u_in,v_in));
    assert(AreCompatible(v_in,w_out));

    u_in.GetDevicePtr()->CrossProduct(u_in.begin(), v_in.begin(), 
        u_in.GetFunLength(), u_in.GetNumVecs(), w_out.begin());
}

template<typename T, enum DeviceType DT>
inline void uyInv(const Vectors<T,DT> &u_in, 
                  const Scalars<T,DT> &y_in, 
                        Vectors<T,DT> &uyInv_out)
{
    assert(AreCompatible(u_in,y_in));
    assert(AreCompatible(y_in,uyInv_out));

    u_in.GetDevicePtr()->uyInv(u_in.begin(), y_in.begin(), 
        u_in.GetFunLength(), u_in.GetNumVecs(), uyInv_out.begin());
}
    
template<typename T, enum DeviceType DT>
inline void avpw(const Scalars<T,DT> &a_in, 
                 const Vectors<T,DT> &v_in, 
                 const Vectors<T,DT> &w_in, 
                       Vectors<T,DT> &avpw_out)
{
    assert(a_in.GetNumFuns() == v_in.GetNumVecs());
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,avpw_out));

    u_in.GetDevicePtr()->avpw(a_in.begin(), v_in.begin(), 
        w_in.begin(), v_in.GetFunLength(), u_in.GetNumVecs(), avpw_out.begin());
}
 
template<typename T, enum DeviceType DT>
inline void xvpw(const Scalars<T,DT> &x_in, 
                 const Vectors<T,DT> &v_in,
                 const Vectors<T,DT> &w_in, 
                       Vectors<T,DT> &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,w_in));
    assert(AreCompatible(w_in,xvpw_out));

    x_in.GetDevicePtr()->xvpw(x_in.begin(), v_in.begin(), 
        w_in.begin(), v_in.GetFunLength(), v_in.GetNumVecs(), xvpw_out.begin());
}

template<typename T, enum DeviceType DT>
inline void xv(const Scalars<T,DT> &x_in, 
               const Vectors<T,DT> &v_in,
                     Vectors<T,DT> &xvpw_out)
{
    assert(AreCompatible(x_in,v_in));
    assert(AreCompatible(v_in,xvpw_out));

    x_in.GetDevicePtr()->xvpw<T>(x_in.begin(), v_in.begin(), 
        NULL, v_in.GetFunLength(), v_in.GetNumVecs(), xvpw_out.begin());
}
    
template<typename T, enum DeviceType DT>
inline void ShufflePoints(const Vectors<T,DT> &u_in, 
                                Vectors<T,DT> &u_out)
{
    assert(AreCompatible(u_in,u_out));
    size_t stride = u_in.GetFunLength();

    size_t dim1 = (u_in.GetOrder() == AxisMajor) ? DIM : stride;
    size_t dim2 = (u_in.GetOrder() == AxisMajor) ? stride : DIM;

    for(size_t ss=0; ss<u_in.GetNumVecs(); ++ss)
        u_in.GetDevicePtr()->Transpose(u_in.begin() + DIM * stride *ss
            , dim1, dim2, u_out.begin() + DIM * stride * ss);
        
    u_out.SetOrder((u_in.GetOrder() == AxisMajor) ? PointMajor : AxisMajor);
}

template<typename T, enum DeviceType DT>
std::ostream& operator<<(std::ostream& output, Scalars<T,DT> &sc)
{
    output<<"=====================================================\n"
          <<"SH order            : "<<sc.GetShOrder()<<"\n"
          <<"Grid size           : ( "
          <<sc.GetGridDim().first<<" , "<<sc.GetGridDim().second<<" )"<<"\n"
          <<"stride              : "<<sc.GetFunLength()<<"\n"
          <<"Number of functions : "<<sc.GetNumFuns()<<"\n"
          <<"=====================================================\n";
    for(typename Scalars<T,DT>::iterator it = sc.begin(); it !=sc.end(); ++it)
    {
        output<<*it<<" ";

        if((it-sc.begin() + 1)%sc.GetGridDim().second == 0)
            output<<endl;
        
        if((it-sc.begin() + 1)%sc.GetFunLength() == 0)
            output<<endl;

    }
    return(output);
}
