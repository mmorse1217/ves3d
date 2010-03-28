#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <cublas.h>
//#include <cutil.h> ///@todo I was not able to find this on Ion.

#define scalar float

extern void cu_trans(scalar *out, scalar *in, int width, int height);

class CudaSht {
private:
  int p;
  int dft_size;
  int leg_mat_size;
  int vesicle_size;

  scalar *leg_temp;
  scalar *trans_in;
  scalar *trans_out;
  scalar *dft_temp;

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

  void post_legendre(scalar *in, scalar *out);
  void pre_legendre(scalar *in, scalar *out);

  void read_leg_mat(scalar **gpu_ptr, char *fname);
//  void cublas_alloc_copy(scalar *cpu_ptr, scalar **gpu_ptr, int rows, int cols);

  void leg_transform(scalar *trans_gpu, scalar *inputs_gpu, scalar *outputs_gpu,
                  int m, int n , int k, int mf, int nf, int kf);
  void back(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs, scalar *trans, scalar *dft);

public:
  void forward(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward_du(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward_dv(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward_d2u(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward_d2v(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
  void backward_duv(scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);

CudaSht(int p, int num_vesicles, char *leg_trans_fname, char *leg_trans_inv_fname,
         char *d1_leg_trans_fname, char *d2_leg_trans_fname, scalar *leg_temp,
         scalar *dft_temp);
~CudaSht();

  static void cublas_alloc_copy(scalar *cpu_ptr, scalar **gpu_ptr, int rows, int cols);
  void test(int n_funs);
};
