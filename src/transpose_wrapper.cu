#define BLOCK_VES3D_DIM 16

extern __global__ void transpose(float *out, float *in, int width, int height);

void cu_trans(float *out, float *in, int width, int height) {
  
  dim3 grid(width / BLOCK_VES3D_DIM + 1, height / BLOCK_VES3D_DIM + 1, 1);
  dim3 threads(BLOCK_VES3D_DIM, BLOCK_VES3D_DIM, 1);
  transpose<<<grid, threads>>>(out, in, width, height);
}
