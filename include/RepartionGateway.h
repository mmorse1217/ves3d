#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

template<typename Container>
class RepartionGateway
{
  private:
    typedef typename Container::value_type value_type;
    
  public:
    typedef void(*GlobalRepart_t)(size_t nv, size_t stride, 
        const value_type* x, const value_type* tension, size_t* nvr, 
        value_type** xr, value_type** tensionr, void* user_ptr);
    
    explicit RepartionGateway(GlobalRepart_t fun_ptr = NULL) :
        g_repart_handle_(fun_ptr) {};
    
    void operator()(Container &coord, 
        Container &tension, void* user_ptr) const
    {
        assert(typeid(Container::getDevice()) == typeid(Device<CPU>)); 
        
        if ( g_repart_handle_ == NULL ) 
            return;
        
        size_t nv(coord.getNumSubs());
        size_t stride(coord.getStride());
        value_type* xr;
        value_type* tensionr;
        size_t nvr(0);
        
        g_repart_handle_(nv, stride, coord.begin(), tension.begin(), &nvr
            &xr, &tensionr, user_ptr);
        
        coord.resize(nvr);
        tension.resize(nvr);
        
        Container::getDevice().Memcpy(coord.begin(), xr, nvr * stride * DIM,
            MemcpyHostToDevice);

        Container::getDevice().Memcpy(tension.begin(), tensionr, nvr * stride,
            MemcpyHostToDevice);
        
        delete[] xr;
        delete[] tensionr;
    }

  private:
    GlobalRepart_t g_repart_handle_;
};
#endif // _REPARTIONGATEWAY_H_
