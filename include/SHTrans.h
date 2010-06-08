#ifndef _SHTRANS_H_
#define _SHTRANS_H_

#include "Device.h"
#include "SHTMats.h"

template<typename Container>
class SHTrans 
{
  public:
    SHTrans(int sh_order_in);
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
    SHTMats<float,CPU> mats_;

    static const float alpha_ = (float) 1.0;
    static const float beta_  = (float) 0.0;

    int p;
    int dft_size;
        
    void DLT(float *trans, const float *inputs, float *outputs, int m, int n , 
        int k, int mf, int nf, int kf) const;
    void back(const float *inputs, float *work_arr, int n_funs, float *outputs, 
        float *trans, float *dft) const;

    void ScaleFreq(const float *shc_in, int n_funs, const float* scaling_coeff, float *shc_out) const;

    float* filter_coeff_;

};

#include "SHTrans.cc"
#endif //_SHTRANS_H_
