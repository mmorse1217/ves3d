#include "cuda_sht.h"
#include "blas_sht.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sys/time.h>
#include <assert.h>

int main(int argc, char ** argv)
{
  srand48(time(0));

  // p is the "order"; we have 2*p points on each "parallel" and p+1
  // points on each meridian
  int p=1;
  if (argc>1)
    p = atoi(argv[1]);

  int dft_size=2*p;

  int num_vesicles = 1;
  if (argc>2)
    num_vesicles = atoi(argv[2]);
  printf("Using p=%d num_vesicles=%d\n",p,num_vesicles);

  cuda_sht  my_sht(p,num_vesicles,"../data/legTrans12", "../data/legTransInv12", "../data/d1legTrans12", "../data/d2legTrans12");
  blas_sht  my_blas_sht(p,num_vesicles,"../data/legTrans12", "../data/legTransInv12", "../data/d1legTrans12", "../data/d2legTrans12");

  int num_dft_inputs = num_vesicles*(p+1);

  // we assume that initially we are given matrix of values in
  // physical coordinates for each vesicle; we assume the matrices are
  // (2*p by p+1) and stored column-wise and contiguosly (matrix after
  // matrix).
  scalar * inputs = (scalar *) malloc (2*p*(p+1)*num_vesicles*sizeof(scalar));
  vector<scalar> input2(2*p*(p+1)*num_vesicles);

  for (int i=0; i<dft_size*num_dft_inputs; i++)
    inputs[i] = input2[i]= drand48();

  scalar * legendre_outputs = (scalar *) malloc (p*(p+2)*num_vesicles*sizeof(scalar));
  vector<scalar> outputs2(p*(p+2)*num_vesicles);

  my_sht.forward(inputs, legendre_outputs);
  my_blas_sht.forward(&input2[0], &outputs2[0]);

//   my_sht.forward(inputs, legendre_outputs);
//   my_blas_sht.forward(&input2[0], &outputs2[0]);
// 
//   my_sht.forward(inputs, legendre_outputs);
//   my_blas_sht.forward(&input2[0], &outputs2[0]);
//   my_blas_sht.forward(&input2[0], &outputs2[0]);
//   my_blas_sht.forward(&input2[0], &outputs2[0]);
//   my_blas_sht.forward(&input2[0], &outputs2[0]);
  

  // compare CPU-BLAS and CUDA-BLAS results
  for(int i=0; i<p*(p+2)*num_vesicles; i++)
    assert(fabs(legendre_outputs[i]-outputs2[i])<1e-5);
  
  my_sht.backward(legendre_outputs, inputs);
  my_blas_sht.backward(&outputs2[0], &input2[0]);

  for(int i=0; i<p*(p+2)*num_vesicles; i++)
    assert(fabs(inputs[i]-input2[i])<1e-5);

  my_sht.forward(inputs, legendre_outputs);
  my_sht.backward(legendre_outputs, &input2[0]);

  for (int i=0; i< 2*p*(p+1); i++)
    assert(fabs(inputs[i]-input2[i])<1e-5);

  
     my_sht.backward_du(legendre_outputs, &input2[0]);
    // same for F_u

     my_sht.backward_dv(legendre_outputs, &input2[0]);
    // F_v
    
     my_sht.backward_d2u(legendre_outputs, &input2[0]);
    // same for F_uu

     my_sht.backward_d2v(legendre_outputs, &input2[0]);
    // same for F_vv
    
     my_sht.backward_duv(legendre_outputs, &input2[0]);
    // same for F_uv


  free(inputs);
  free(legendre_outputs);
  return 0;
}
