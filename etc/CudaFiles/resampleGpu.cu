#include <cuda.h>

void ResampleGpu(int p, int n_funs, int q, float *shc_p, float *shc_q) {

  float *out_deb = shc_q;
  int leg_order = p + 1;
  int new_leg_order = q + 1;
  int min_leg_order = (leg_order < new_leg_order) ? leg_order : new_leg_order;

  for(int v = 0; v < n_funs; v++) {
    cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
    shc_q += min_leg_order;
    shc_p += min_leg_order;
    if (new_leg_order > leg_order) {
      cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
      shc_q += (new_leg_order - leg_order);
    }
    if (leg_order > new_leg_order)
      shc_p += (leg_order - new_leg_order);
  }
  leg_order--;
  new_leg_order--;
  min_leg_order--;

  for(; min_leg_order > 1; min_leg_order--, leg_order--, new_leg_order--) {
    for(int v = 0; v < n_funs; v++) {
      cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
      shc_q += min_leg_order;
      shc_p += min_leg_order;
      if (new_leg_order > leg_order) {
        cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
        shc_q += (new_leg_order - leg_order);
      }
      if (leg_order > new_leg_order)
        shc_p += (leg_order - new_leg_order);
    }
    for(int v = 0; v < n_funs; v++) {
      cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
      shc_q += min_leg_order;
      shc_p += min_leg_order;
      if (new_leg_order > leg_order) {
        cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
        shc_q += (new_leg_order - leg_order);
      }
      if (leg_order > new_leg_order)
        shc_p += (leg_order - new_leg_order);
    }
  }

  for(int v = 0; v < n_funs; v++) {
    cudaMemcpy(shc_q, shc_p, sizeof(float) * min_leg_order, cudaMemcpyDeviceToDevice);
    shc_q += min_leg_order;
    shc_p += min_leg_order;
    if (new_leg_order > leg_order) {
      cudaMemset(shc_q, 0, sizeof(float) * (new_leg_order - leg_order));
      shc_q += (new_leg_order - leg_order);
    }
    if (leg_order > new_leg_order)
      shc_p += (leg_order - new_leg_order);
  }

  leg_order--;
  new_leg_order--;
  min_leg_order--;

  float *outputs_end = out_deb + n_funs * q * (q + 2);
  if (shc_q < outputs_end) {
    cudaMemset(shc_q, 0, sizeof(float) * (outputs_end - shc_q));
  }
}
