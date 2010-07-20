#ifndef _REPARTIONGATEWAY_H_
#define _REPARTIONGATEWAY_H_

template<typename Container>
class RepartionGateway
{
  private:
    typedef typename Container::value_type value_type;
    
  public:
    typedef void(*GlobalRepart_t)(size_t np, value_type* x, value_type* tension,
        size_t* npr, value_type** xr, value_type** tensionr, void*);
    
    RepartitionGateway(GlobalRepart_t fun_ptr);
    
    void operator()(Container &coord, Container &tension, void* user_ptr) const
    {
        g_repart_handle_(...);
    }

  private:
    GlobalRepart_t g_repart_handle_;
};
#endif // _REPARTIONGATEWAY_H_
