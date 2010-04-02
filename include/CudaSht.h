#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cublas.h>
#include "transpose_kernel.h"
//#include <cutil.h> ///@todo I was not able to find this on Ion.
#define scalar float

class CudaSht {
  public:
    int p;
    int dft_size;
    int leg_mat_size;
    int vesicle_size;

    scalar *leg_trans;
    scalar *leg_trans_inv;
    scalar *d1_leg_trans;
    scalar *d2_leg_trans;
    scalar *dft_forward;
    scalar *dft_backward;
    scalar *dft_d1backward;
    scalar *dft_d2backward;

    void gen_dft_forward();
    void gen_dft_backward();
    void gen_dft_d1backward();
    void gen_dft_d2backward();

    void post_legendre(const scalar *in, scalar *out);
    void pre_legendre(const scalar *in, scalar *out);

    void leg_transform(scalar *trans_gpu, const scalar *inputs_gpu, scalar *outputs_gpu,
        int m, int n , int k, int mf, int nf, int kf);
    void back(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs, scalar *trans, scalar *dft);

  public:
    CudaSht();
    CudaSht(int p, scalar *dft_forward, scalar *dft_backward, 
        scalar *dft_d1backward, scalar *dft_d2backward,
        scalar *leg_trans, scalar *leg_trans_inv, scalar* d1_leg_trans,
        scalar *d2_leg_trans);

    ~CudaSht();
    
    void InitializeCudaSht(int p, int num_vesicles, char *leg_trans_fname, 
        char *leg_trans_inv_fname, char *d1_leg_trans_fname, 
        char *d2_leg_trans_fname);

    void forward(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_du(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_dv(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_d2u(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_d2v(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_duv(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);

    void InitializeCudaSht(int p, scalar *dft_forward, scalar *dft_backward, 
        scalar *dft_d1backward, scalar *dft_d2backward,
        scalar *leg_trans, scalar *leg_trans_inv, scalar* d1_leg_trans,
        scalar *d2_leg_trans);

    static void cublas_alloc_copy(scalar *cpu_ptr, scalar **gpu_ptr, int rows, int cols);
    void test(int n_funs);
};
