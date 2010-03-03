/**
 * @file   SHTrans.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Jan 31 13:50:30 2010
 * 
 * @brief  The implementation of the SHTrans class. 
 */

template<typename T> 
SHTrans<T>::SHTrans(Device<T> &device_in) : 
    device_(device_in), p_(0), n_funs_(0), shc_(0){}

template<typename T> 
SHTrans<T>::SHTrans(Device<T> &device_in, int p_in, int n_funs_) : 
    device_(device_in), p_(p_in),   n_funs_(n_funs_),
    cudaTransClass(p_,n_funs_,"../data/legTrans12", 
        "../data/legTransInv12", "../data/d1legTrans12",
        "../data/d2legTrans12")
{
    int size = p_*(p_+2)*n_funs_;
    shc_ = device_.Malloc(size);
}

template<typename T> 
SHTrans<T>::~SHTrans()
{
    device_.Free(shc_);
    shc_ = 0;
}


template <typename T> 
void SHTrans<T>::AllDerivatives(const Scalars<T> &f_in, 
    Scalars<T> &Duf_out, Scalars<T> &Dvf_out, 
    Scalars<T> &Duuf_out, Scalars<T> &Duvf_out, 
    Scalars<T> &Dvvf_out)
{
    assert(    f_in.device_ == device_);
    assert( Duf_out.device_ == device_);
    assert( Dvf_out.device_ == device_);
    assert(Duuf_out.device_ == device_);
    assert(Duvf_out.device_ == device_);
    assert(Dvvf_out.device_ == device_);

    cudaTransClass.forward(f_in.data_, shc_);
    cudaTransClass.backward_du (shc_,  Duf_out.data_);
    cudaTransClass.backward_dv (shc_,  Dvf_out.data_);
    cudaTransClass.backward_d2u(shc_, Duuf_out.data_);
    cudaTransClass.backward_duv(shc_, Duvf_out.data_);
    cudaTransClass.backward_d2v(shc_, Dvvf_out.data_);
}

template <typename T> 
void SHTrans<T>::FirstDerivatives(const Scalars<T> &f_in, 
    Scalars<T> &Duf_out, Scalars<T> &Dvf_out) 
{
    assert(    f_in.device_ == device_);
    assert( Duf_out.device_ == device_);
    assert( Dvf_out.device_ == device_);
    
    cudaTransClass.forward(f_in.data_, shc_);
    cudaTransClass.backward_du(shc_, Duf_out.data_);
    cudaTransClass.backward_dv(shc_, Dvf_out.data_);
}




