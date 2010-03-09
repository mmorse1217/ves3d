#include <fstream>
#include <iostream>
#include <math.h>
#include <mkl.h>
#include <stdlib.h>

#define scalar float

class BlasSht {
  public:
    int p;
    //int num_vesicles;
    //int num_dft_inputs;
    int dft_size;
    int leg_mat_size;
    int vesicle_size;
    scalar alpha;
    scalar beta;
    
    //scalar *trans_in;
    //scalar *trans_out;
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

    void read_leg_mat(scalar *leg_ptr, char *fname);

    void leg_transform(scalar *trans, const scalar *inputs, scalar *outputs,
        int m, int n , int k, int mf, int nf, int kf);
    void back(const scalar *inputs, scalar *work_arr, int n_funs,
               scalar *outputs, scalar *trans, scalar *dft);
    void transpose(scalar *out, scalar *in, int width, int height);

  public:
    BlasSht();

    BlasSht(int p, char *leg_trans_fname, char *leg_trans_inv_fname,
        char *d1_leg_trans_fname, char *d2_leg_trans_fname, scalar *dft_forward,
        scalar *dft_backward, scalar *dft_d1backward, scalar *dft_d2backward,
        scalar *leg_trans, scalar *leg_trans_inv, scalar* d1_leg_trans,
        scalar *d2_leg_trans);
    ~BlasSht();

    void forward(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_du(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_dv(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_d2u(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_d2v(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);
    void backward_duv(const scalar *inputs, scalar *work_arr, int n_funs, scalar *outputs);

    void InitializeBlasSht(int p, char *leg_trans_fname,
        char *leg_trans_inv_fname, char *d1_leg_trans_fname, 
        char *d2_leg_trans_fname, scalar *dft_forward, scalar *dft_backward, 
        scalar *dft_d1backward, scalar *dft_d2backward, scalar *leg_trans, 
        scalar *leg_trans_inv, scalar* d1_leg_trans, scalar *d2_leg_trans);

    void test(int n_funs);
};
