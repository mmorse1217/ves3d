#include <cuda.h>

void CircShiftGpu(const float *arr_in, int n_vecs, int vec_length, int shift, float *arr_out) {
  shift = shift % vec_length;
  if (shift < 0) {
    shift += vec_length;
  }
  int base_in, base_out;
  for (int ii = 0; ii < n_vecs; ii++) {
    base_out = ii * vec_length;
    base_in = base_out + vec_length - shift;
    cudaMemcpy(arr_out + base_out, arr_in + base_in, sizeof(float) * shift, cudaMemcpyDeviceToDevice);
    base_in = base_out;
    base_out += shift;
    cudaMemcpy(arr_out + base_out, arr_in + base_in, sizeof(float) * (vec_length - shift),
               cudaMemcpyDeviceToDevice);
  }
}
