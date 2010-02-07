/**
 * @file   SphHarm.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Sun Jan 31 13:50:30 2010
 * 
 * @brief  The implementation of the SphHarm class. 
 */

template<typename scalarType> 
SphHarm<scalarType>::SphHarm() : 
  p_(0), number_of_functions_(0), shc_(NULL){}


template<typename scalarType> 
SphHarm<scalarType>::SphHarm(int p_in, int number_of_functions_) : 
  p_(p_in),   number_of_functions_(number_of_functions_),
  cudaTransClass(p_,number_of_functions_,"../data/legTrans12", 
		 "../data/legTransInv12", "../data/d1legTrans12",
		 "../data/d2legTrans12")
{
  int size = p_*(p_+2)*number_of_functions_;
  shc_ = new scalarType[size];
}

template<typename scalarType> 
SphHarm<scalarType>::~SphHarm()
{
  delete[] shc_;
  shc_ = NULL;
}


template <typename scalarType> 
void SphHarm<scalarType>::Derivatives(const SHScalars<scalarType>& f_in, 
    SHScalars<scalarType>& Duf_out, SHScalars<scalarType>& Dvf_out, 
    SHScalars<scalarType>& Duuf_out, SHScalars<scalarType>& Duvf_out, 
    SHScalars<scalarType>& Dvvf_out)
{
  cudaTransClass.forward(f_in.data_, shc_);
  cudaTransClass.backward_du (shc_, Duf_out.data_ );
  cudaTransClass.backward_dv (shc_, Dvf_out.data_ );
  cudaTransClass.backward_d2u(shc_, Duuf_out.data_);
  cudaTransClass.backward_duv(shc_, Duvf_out.data_);
  cudaTransClass.backward_d2v(shc_, Dvvf_out.data_);
}

