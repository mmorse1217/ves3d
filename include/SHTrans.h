#ifndef _SHTRANS_H_
#define _SHTRANS_H_

#include "Device.h"
#include "SHTMats.h"

template<typename T, enum DeviceType DT>
class SHTrans {
  private:
    Device<DT> *device_;
    SHTMats<T, DT> mats_;

    T alpha_;
    T beta_;

    int p;
    int dft_size;
    int leg_mat_size;
        
    void DLT(T *trans, const T *inputs, T *outputs, int m, int n , int k, int mf, int nf, int kf);
    void back(const T *inputs, T *work_arr, int n_funs, T *outputs, T *trans, T *dft);
    
  public:
    T *leg_trans;
    T *leg_trans_inv;
    T *d1_leg_trans;
    T *d2_leg_trans;
    T *dft_forward;
    T *dft_backward;
    T *dft_d1backward;
    T *dft_d2backward;

    SHTrans(Device<DT> *dev, int p_in);
    ~SHTrans();

    void forward(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward_du(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward_dv(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward_d2u(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward_d2v(const T *inputs, T *work_arr, int n_funs, T *outputs);
    void backward_duv(const T *inputs, T *work_arr, int n_funs, T *outputs);
};

#include "SHTrans.cc"
#endif //_SHTRANS_H_
