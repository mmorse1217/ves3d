#ifndef _SHTRANS_H_
#define _SHTRANS_H_

#include "Device.h"
#include "SHTMats.h"
#include "OperatorsMats.h"

template<typename Container>
class SHTrans 
{
  private:
    typedef typename Container::value_type value_type;

  public:
    SHTrans(int sh_order_in, OperatorsMats<value_type> &mats, int filter_freq = -1);
    ~SHTrans();

    void FirstDerivatives(const Container &in, Container &work, Container &shc, Container &du, Container &dv) const;

    void forward(const Container &in, Container &work, Container &shc) const;
   
    void backward    (const Container &shc, Container &work, Container &out) const;
    void backward_du (const Container &shc, Container &work, Container &out) const;
    void backward_dv (const Container &shc, Container &work, Container &out) const;
    void backward_d2u(const Container &shc, Container &work, Container &out) const;
    void backward_d2v(const Container &shc, Container &work, Container &out) const;
    void backward_duv(const Container &shc, Container &work, Container &out) const;

    void Filter(const Container &in, Container &work, Container &shc, Container &out) const;
    
  private:
    const Device<CPU> *device_;
    SHTMats<value_type,CPU> mats_;
    
    static const value_type alpha_ = (value_type) 1.0;
    static const value_type beta_  = (value_type) 0.0;

    int p;
    int dft_size;
        
    void DLT(value_type *trans, const value_type *inputs, value_type *outputs, int m, int n , 
        int k, int mf, int nf, int kf) const;
    void back(const value_type *inputs, value_type *work_arr, int n_funs, value_type *outputs, 
        value_type *trans, value_type *dft) const;

    void ScaleFreq(const value_type *shc_in, int n_funs, const value_type* scaling_coeff, value_type *shc_out) const;

    value_type* filter_coeff_;
};

#include "SHTrans.cc"
#endif //_SHTRANS_H_
