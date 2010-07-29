#include <iostream>

#include <cublas.h>
#include <math.h>
#include <stdlib.h>

void tfqmr(float *x, int n, float *b, float rres, int kmax, float *scratch,
            void (*atv)(int n, float*, float*)) {

  float *y[2];
  float *u[2];
  float *ax, *w, *d, *v, *r, *t;

  y[0] = scratch;
  y[1] = y[0] + n;
  u[0] = y[1] + n;
  u[1] = u[0] + n;
  ax = u[1] + n;
  w = ax + n;
  d = w + n;
  v = d + n;
  r = v + n;
  t = r + n;

  float errtol = rres * cublasSnrm2(n, b, 1);

  atv(n, x, ax);

  cublasScopy(n, b, 1, r, 1);
  cublasSaxpy(n, -1.0F, ax, 1, r, 1);

  cublasScopy(n, r, 1, w, 1);
  cublasScopy(n, r, 1, y[0], 1);
  atv(n, y[0], v);
  cublasScopy(n, v, 1, u[0], 1);
  float theta = 0;
  float eta = 0;
  float tau = cublasSnrm2(n, r, 1);
  float rho = tau * tau;

  int k = 0;
  while (k < kmax) {
    k++;
    float sigma = cublasSdot(n, r, 1, v, 1);
    if(sigma == 0.0F) {
      std::cout << "TFQMR breakdown, sigma=0";
      abort();
    }
    float alpha = rho / sigma;

    for (int j = 0; j < 2; j++) {
      if (j == 1) {
        cublasScopy(n, y[0], 1, y[1], 1);
        cublasSaxpy(n, -alpha, v, 1, y[1], 1);
        atv(n, y[1], u[1]);
      }
      int m = 2 * (k + 1) - 2 + (j + 1);
      cublasSaxpy(n, -alpha, u[j], 1, w, 1);
      cublasScopy(n, y[j], 1, t, 1);
      cublasSaxpy(n, theta * theta * eta / alpha, d, 1, t, 1);
      cublasScopy(n, t, 1, d, 1);

      theta = cublasSnrm2(n, w, 1) / tau;
      float c = 1 / sqrtf(1 + theta * theta);
      tau = tau * theta * c;
      eta = c * c * alpha;
      cublasSaxpy(n, eta, d, 1, x, 1);

      if (tau * sqrtf(m +1) <= errtol) {
        return;
      }
    }
    if (rho == 0.0F) {
      std::cout<< "TFQMR breakdown, rho=0";
      abort();
    }

    float rhon = cublasSdot(n, r, 1, w, 1);
    float beta = rhon / rho;
    rho = rhon;
    cublasScopy(n, w, 1, y[0], 1);
    cublasSaxpy(n, beta, y[1], 1, y[0], 1);
    atv(n, y[0], u[0]);
    cublasSaxpy(n, beta, v, 1, u[1], 1);
    cublasScopy(n, u[0], 1, v, 1);
    cublasSaxpy(n, beta, u[1], 1, v, 1);
  }
}
