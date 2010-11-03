#define BLOCK_DIM 16

extern __global__ void transpose(float *out, float *in, int width, int height);

void cu_trans(float *out, float *in, int width, int height) {
  
  dim3 grid(width / BLOCK_DIM + 1, height / BLOCK_DIM + 1, 1);
  dim3 threads(BLOCK_DIM, BLOCK_DIM, 1);
  transpose<<<grid, threads>>>(out, in, width, height);
}
