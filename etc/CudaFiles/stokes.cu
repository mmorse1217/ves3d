#include <cuda.h>
#include <stdio.h>

#define BLOCK_HEIGHT 128
#define scalar float
#define scalar3 float3
#define PI_8I 0.0397887358F


__global__
void stokes(int m, int n, int t_head, scalar *T, scalar *S, scalar *D, scalar *U, scalar *Q) {
  scalar3 trg_reg;
  scalar3 src_reg;
  scalar3 pot_reg;
  scalar3 dis_reg;
  scalar3 u_reg;
  __shared__
  scalar3 u_sh[BLOCK_HEIGHT];

  int t_off = blockIdx.x * 3 * m + t_head + blockIdx.y;

  trg_reg.x = T[t_off];
  trg_reg.y = T[m + t_off];
  trg_reg.z = T[m + m + t_off];

  u_reg = make_float3(0.0, 0.0, 0.0);

  int block_off = blockIdx.x * 3 * m + threadIdx.x;
  int s = threadIdx.x;

  while(block_off < blockIdx.x * 3 * m + m) {
    src_reg.x = S[block_off];
    src_reg.y = S[block_off + m];
    src_reg.z = S[block_off + m + m];

    pot_reg.x = D[block_off] * Q[s];
    pot_reg.y = D[block_off + m] * Q[s];
    pot_reg.z = D[block_off + m + m] * Q[s];

    dis_reg.x = src_reg.x - trg_reg.x;
    dis_reg.y = src_reg.y - trg_reg.y;
    dis_reg.z = src_reg.z - trg_reg.z;

    scalar inv_r = rsqrtf(dis_reg.x * dis_reg.x + dis_reg.y * dis_reg.y
                          + dis_reg.z * dis_reg.z);

    scalar tmp_scal = (dis_reg.x * pot_reg.x + dis_reg.y * pot_reg.y
                       + dis_reg.z * pot_reg.z) * inv_r * inv_r;
    pot_reg.x += tmp_scal * dis_reg.x;
    pot_reg.y += tmp_scal * dis_reg.y;
    pot_reg.z += tmp_scal * dis_reg.z;

    u_reg.x += pot_reg.x * inv_r;
    u_reg.y += pot_reg.y * inv_r;
    u_reg.z += pot_reg.z * inv_r;

    block_off += BLOCK_HEIGHT;
    s += BLOCK_HEIGHT;
  }

  u_sh[threadIdx.x].x = u_reg.x;
  u_sh[threadIdx.x].y = u_reg.y;
  u_sh[threadIdx.x].z = u_reg.z;

  int off = 1;
  int stride = 2;
  while (off != BLOCK_HEIGHT) {
    if (threadIdx.x % stride == 0) {
      syncthreads();
      u_sh[threadIdx.x].x += u_sh[threadIdx.x + off].x;
      u_sh[threadIdx.x].y += u_sh[threadIdx.x + off].y;
      u_sh[threadIdx.x].z += u_sh[threadIdx.x + off].z;
    }
    off = stride;
    stride *= 2;
  }
  if (threadIdx.x == 0) {
    U[t_off] = u_sh[0].x * PI_8I;
    U[m + t_off] = u_sh[0].y * PI_8I;
    U[m + m + t_off] = u_sh[0].z * PI_8I;
  }

}


void cuda_stokes(int m, int n, int t_head, int t_tail, scalar *T, scalar *S, scalar *D, scalar *U, scalar *Q) {
  dim3 grid;
  grid.x = n;
  grid.y = t_tail - t_head;
  stokes<<<grid, BLOCK_HEIGHT>>> (m, n, t_head, T, S, D, U, Q);
}
