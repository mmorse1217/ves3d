#include <cuda.h>

#define BLOCK_HEIGHT 128

__global__
void xyMKernel(const float *x_in, const float *y_in, int length, int stride, float *xy_out) {
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < length) {
    xy_out[idx] = x_in[idx] * y_in[idx % stride];
  }
}

void xyMGpu(const float* x_in, const float *y_in, int stride, int num_surfs, float *xy_out) {
  int length = stride * num_surfs;
  int grid = length / BLOCK_HEIGHT + 1;
  xyMKernel<<<grid, BLOCK_HEIGHT>>> (x_in, y_in, length, stride, xy_out);
}

void ScaleFreqsGpu(int p, int n_funs, float *shc_in, float *alpha, float *shc_out) {
    int leg_order = p+1;
    xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
    alpha += leg_order;
    shc_in += n_funs * leg_order;
    shc_out += n_funs * leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
        alpha += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
        
        // then process sine
        xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
        alpha += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
    }
    
    // process last cosine
    xyMGpu(shc_in, alpha, leg_order, n_funs, shc_out);
}
