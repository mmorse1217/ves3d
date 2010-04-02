#include "CudaSht.h"
#include "DataIO.h"
#include "OperatorsMats.h"
#include "cuda_runtime.h"
//#include "BlasSht.h"

#include "DeviceCPU.h"

using namespace std;
typedef float T;

int main(int argc, char *argv[])
{
    int p = 12;
    int leg_size = (p + 1) * (p + 1) * (p + 2);

    CudaSht shtCuda;
    BlasSht shtBlas;

    //reading mats
    DeviceCPU<float> gpu_f;
    DataIO<float> gpuIO(gpu_f," ", 0);

    bool readFromFile = true;
    OperatorsMats<T> mats(gpuIO,p, 2*p, readFromFile);
    
    //Initializing Blas
    T *dft_forward_host    = (T*) malloc(4 * p * p * sizeof(T));
    T *dft_backward_host   = (T*) malloc(4 * p * p * sizeof(T));
    T *dft_d1backward_host = (T*) malloc(4 * p * p * sizeof(T));
    T *dft_d2backward_host = (T*) malloc(4 * p * p * sizeof(T));

    shtBlas.InitializeBlasSht(p, dft_forward_host, dft_backward_host, 
        dft_d1backward_host, dft_d2backward_host, mats.leg_trans_p_, 
        mats.leg_trans_inv_p_, mats.d1_leg_trans_p_, mats.d2_leg_trans_p_);

    //Initializing Cuda
    T *dft_forward;    
    T *dft_backward;
    T *dft_d1backward;
    T *dft_d2backward;

    T *leg_trans;
    T *leg_trans_inv;
    T *d1_leg_trans;
    T *d2_leg_trans;

    cudaMalloc(&dft_forward   , 4 * p * p * sizeof(T));
    cudaMalloc(&dft_backward  , 4 * p * p * sizeof(T));
    cudaMalloc(&dft_d1backward, 4 * p * p * sizeof(T));
    cudaMalloc(&dft_d2backward, 4 * p * p * sizeof(T));

    cudaMalloc(&leg_trans    , leg_size * sizeof(T));
    cudaMalloc(&leg_trans_inv, leg_size * sizeof(T));
    cudaMalloc(&d1_leg_trans , leg_size * sizeof(T));
    cudaMalloc(&d2_leg_trans , leg_size * sizeof(T));
    
    cudaMemcpy(leg_trans    , mats.leg_trans_p_    , leg_size* sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(leg_trans_inv, mats.leg_trans_inv_p_, leg_size* sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d1_leg_trans , mats.d1_leg_trans_p_ , leg_size* sizeof(T), cudaMemcpyHostToDevice);
    cudaMemcpy(d2_leg_trans , mats.d2_leg_trans_p_ , leg_size* sizeof(T), cudaMemcpyHostToDevice);

    shtCuda.InitializeCudaSht(p, dft_forward, dft_backward, dft_d1backward, dft_d2backward,
        leg_trans, leg_trans_inv, d1_leg_trans, d2_leg_trans);

    //Testing
    int n_funs = 1;
    int length = 2 * p * (p + 1) * n_funs;
    int sh_len = p * (p + 2) * n_funs;

    T *x_ref = (T*) malloc(length * sizeof(T));
    T *work_arr_ref = (T*) malloc( 2 * length * sizeof(T));
    T *shc_ref = (T*) malloc( sh_len * sizeof(T));

    srand48(time(0));
    for(int idx=0;idx<length;idx++)
            x_ref[idx] = (T) drand48();

    shtBlas.forward(x_ref, work_arr_ref, n_funs, shc_ref);
    
    T *x;
    T *work_arr;
    T *shc;
    T *shc_host = (T*) malloc( sh_len * sizeof(T));

    cudaMalloc(&x , length * sizeof(T));
    cudaMalloc(&work_arr , 2 * length * sizeof(T));
    cudaMalloc(&shc , sh_len * sizeof(T));

    cudaMemset(shc, 1, sizeof(float)* sh_len);
    cudaMemcpy(x , x_ref , length* sizeof(T), cudaMemcpyHostToDevice);
    shtCuda.forward(x, work_arr, n_funs, shc);
    cudaMemcpy(shc_host , shc , sh_len* sizeof(T), cudaMemcpyDeviceToHost);

    for(int idx=0;idx<sh_len;idx++)
        cout<<shc_ref[idx]-shc_host[idx]<<endl;

    //Freeing
    free(dft_forward_host);    
    free(dft_backward_host);
    free(dft_d1backward_host);
    free(dft_d2backward_host);
    
    cudaFree(dft_forward);    
    cudaFree(dft_backward);
    cudaFree(dft_d1backward);
    cudaFree(dft_d2backward);

    cudaFree(leg_trans);
    cudaFree(leg_trans_inv);
    cudaFree(d1_leg_trans);
    cudaFree(d2_leg_trans);

    free(x_ref);
    free(work_arr_ref);
    free(shc_ref);
    
    cudaFree(x);
    cudaFree(work_arr);
    cudaFree(shc);
    free(shc_host);

    return 0;
}
