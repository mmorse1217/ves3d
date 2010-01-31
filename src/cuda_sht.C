#include "cuda_sht.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sys/time.h>
#include <assert.h>
#include "cublas.h"

#include <iostream>
#include <fstream>
#include "sstream"
#include <vector>

using namespace std;

// constructor
cuda_sht::cuda_sht(int pp, int nnum_vesicles, char * legTransFname, char* legTransInvFname, char* d1legTransFname, char* d2legTransFname)
{
  p=pp;
  num_vesicles = nnum_vesicles;
  int dft_size = 2*p;

  /* ====== load Legendre matrices,  create DFT matrices ======  */
  // load forward Legendre matrices
  int leg_mat_size = (p+1)*(p+1)*(p+2);
  leg_mat_size = leg_mat_size/2;
  legTrans.resize(leg_mat_size);
  ifstream leg_forward_file (legTransFname);
  if (leg_forward_file.is_open())
  {
    int idx=0;
    while (idx<leg_mat_size)
    {
      leg_forward_file>>legTrans[idx++];
      // cout<<legendre_matrices[idx-1]<<endl;
    }
    leg_forward_file.close();
  }
  else 
  {
    cout << "Unable to open file"; 
    abort();
  }
  
  // DFT matrix for analysis
  dft_forward.resize(dft_size*dft_size);
  for(int j=0; j<dft_size; j++)
    dft_forward[dft_size*j] = 1.0/dft_size;
  
  for(int j=0; j<dft_size; j++) 
    for(int i=1; i<p; i++) 
    {
      dft_forward[2*i-1+dft_size*j] = cos(M_PI*i*j/p)/dft_size*2;  // 2*pi*i*j/(2*p)
      dft_forward[2*i+dft_size*j] = sin(M_PI*i*j/p)/dft_size*2;  // 2*pi*i*j/(2*p)
    }

  // last cosine (sawtooth)
  for(int j=0; j<dft_size; j++)
    dft_forward[dft_size-1 + dft_size*j] = cos(M_PI*j)/dft_size;

  // DFT matrix for synthesis  (will be used transposed)
  dft_backw.resize(2*p*2*p);

  for(int j=0; j<dft_size; j++)
    dft_backw[dft_size*j] = 1.0;
  
  for(int j=0; j<dft_size; j++) 
    for(int i=1; i<p; i++) 
    {
      dft_backw[2*i-1+dft_size*j] = cos(M_PI*i*j/p);
      dft_backw[2*i+dft_size*j] = sin(M_PI*i*j/p);
    }

  // last cosine (sawtooth)
  for(int j=0; j<dft_size; j++)
    dft_backw[dft_size-1 + dft_size*j] = cos(M_PI*j);

  // load backward Legendre matrices
  legTransInv.resize(leg_mat_size);
  ifstream leg_backw_file (legTransInvFname);
  if (leg_backw_file.is_open())
  {
    int idx=0;
    while (idx<leg_mat_size)
      leg_backw_file>>legTransInv[idx++];
    leg_backw_file.close();
  }
  else 
  {
    cout << "Unable to open file\n"; 
    abort();
  }

  // load first derivatives of backward Legendre
  d1legTrans.resize(leg_mat_size);
  ifstream leg_d1backw_file (d1legTransFname);
  if (leg_d1backw_file.is_open())
  {
    int idx=0;
    while (idx<leg_mat_size)
      leg_d1backw_file>>d1legTrans[idx++];
    leg_d1backw_file.close();
  }
  else 
  {
    cout << "Unable to open file\n"; 
    abort();
  }

  // load second derivatives of backward Legendre
  d2legTrans.resize(leg_mat_size);
  ifstream leg_d2backw_file (d2legTransFname);
  if (leg_d2backw_file.is_open())
  {
    int idx=0;
    while (idx<leg_mat_size)
      leg_d2backw_file>>d2legTrans[idx++];
    leg_d2backw_file.close();
  }
  else 
  {
    cout << "Unable to open file\n"; 
    abort();
  }

  dft_d1backw.resize(4*p*p); // first derivative of Fourier basis; to be used transposed for synthesis

  for(int j=0; j<dft_size; j++)
    dft_d1backw[dft_size*j] = 0;
  
  for(int j=0; j<dft_size; j++) 
    for(int i=1; i<p; i++) 
    {
      dft_d1backw[2*i-1+dft_size*j] = -i*sin(M_PI*i*j/p);
      dft_d1backw[2*i+dft_size*j] =    i*cos(M_PI*i*j/p);
    }

  // last cosine (sawtooth) turns into sine which is zero at all points
  for(int j=0; j<dft_size; j++)
    dft_d1backw[dft_size-1 + dft_size*j] = 0;


  dft_d2backw.resize(4*p*p); // second derivative of Fourier basis; to be used transposed for synthesis

  for(int j=0; j<dft_size; j++)
    dft_d2backw[dft_size*j] = 0;
  
  for(int j=0; j<dft_size; j++) 
    for(int i=1; i<p; i++) 
    {
      dft_d2backw[2*i-1+dft_size*j] = -i*i*cos(M_PI*i*j/p);
      dft_d2backw[2*i+dft_size*j] = -i*i*sin(M_PI*i*j/p);
    }

  // last cosine (sawtooth) 
  for(int j=0; j<dft_size; j++)
    dft_d2backw[dft_size-1 + dft_size*j] = -p*p*cos(M_PI*j);

  // allocate temporary storage for computations
  temp_data.resize(2*p*(p+1)*num_vesicles);
  temp_data2.resize(2*p*(p+1)*num_vesicles);

  // initialize  cudaBlas library
  cublasInit();
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
}

cuda_sht::~cuda_sht()
{
  cublasShutdown();
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
}

void cuda_sht::forward(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  int dft_size=2*p;
  int num_dft_inputs = num_vesicles*(p+1);

  // we assume that initially we are given matrix of values in physical coordinates for each vesicle; we assume the matrices are (2*p by p+1) and stored column-wise and contiguosly (matrix after matrix).

  // First compute dft transforms (along each fixed latitude)

  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
    scalar * devPtrA;
    scalar * devPtrB;
    scalar * devPtrC;
    scalar * dft_outputs;
    cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)&dft_forward[0], dft_size, (void*)devPtrA, dft_size);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), inputs, dft_size, devPtrB, dft_size);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


    // actual computation:
#ifdef CUDA_SHT_PROFILING
    double ss = get_seconds();
#endif
  for (int ii=0; ii<1; ii++)
  {
    cublasSgemm ('N', 'N', dft_size, num_dft_inputs,
	dft_size, 1.0F, devPtrA, dft_size,
	devPtrB, dft_size, 0.0F,
	devPtrC, dft_size);
  }
    ss = get_seconds()-ss ;
    cout<<"CUDA BLAS Forward DFT time: "<<ss<< endl;
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    dft_outputs = &temp_data[0];
    cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


 //  ----- rearrange data for Legendre transform -----
  scalar * legendre_inputs = &temp_data2[0];
  int leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	legendre_inputs[leg_input_pointer++] = dft_outputs[vesicle_size*vesicle_num + dft_size*i + freq];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);
  // ---- end of rearranging inputs ----

  scalar * legendre_outputs =  outputs;
  int leg_outputs_pointer = 0;
  leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = & legTrans[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1-freq)*(p+1), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1-freq, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1-freq, (void*)devPtrA, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1-freq)*(p+1);

    cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1);


    // actual computation:
    cublasSgemm ('N', 'N', p+1-freq, num_legendre_inputs,
	p+1, 1.0F, devPtrA, p+1-freq,
	devPtrB, p+1, 0.0F,
	devPtrC, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


    cublasGetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), devPtrC, p+1-freq, &legendre_outputs[leg_outputs_pointer], p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1-freq)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], 1, (void*)devPtrA, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);
  assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', 1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, 1,
      devPtrB, p+1, 0.0F,
      devPtrC, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (1, num_legendre_inputs, sizeof(scalar), devPtrC, 1, legendre_outputs+leg_outputs_pointer, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += num_legendre_inputs;
  assert(leg_outputs_pointer == p*(p+2)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

};


void cuda_sht::backward(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &legTransInv[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
#ifdef CUDA_SHT_PROFILING
    double ss = get_seconds();
#endif
  for (int ii=0; ii<1; ii++)
  {
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  }
    ss = get_seconds()-ss ;
    cout<<"CUDA BLAS Backward DFT time: "<<ss<< endl;
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

void cuda_sht::backward_du(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &d1legTrans[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

void cuda_sht::backward_d2u(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &d2legTrans[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

void cuda_sht::backward_dv(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &legTransInv[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_d1backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

void cuda_sht::backward_d2v(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &legTransInv[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_d2backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

void cuda_sht::backward_duv(scalar * inputs, scalar * outputs)
{
  // p is the "order"; we have 2*p points on each "parallel" and p+1 points on each meridian

  // we first do backward Legendre transform, then backward DFT
  scalar * legendre_outputs =  &temp_data[0];
  scalar * legendre_inputs = inputs;
  int leg_outputs_pointer = 0;
  int leg_input_pointer = 0;

  int leg_matrix_pointer = 0;

  // we have even-order real dft; this means we don't have first sine (sine of zero frequency) and last sine (sine of half-order frequency)
  // -------- process zeroth frequency (constant) ---------
  int num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  scalar * devPtrA;
  scalar * devPtrB;
  scalar * devPtrC;

  cublasAlloc ((p+1)*(p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * legendre_matrices = &d1legTrans[0];
  cublasSetMatrix (p+1, p+1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer] , p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1)*(p+1);
  
  cublasSetMatrix (p+1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1, devPtrB, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs*(p+1);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      p+1, 1.0F, devPtrA, p+1,
      devPtrB, p+1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // process remaining frequencies except the last cosine
  for (int freq=1; freq<p; freq++)
  {
    num_legendre_inputs = 2*num_vesicles; // sine and cosine for each vesicle

    // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
    cublasAlloc ((p+1)*(p+1-freq), sizeof(scalar), (void**)&devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1-freq)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasSetMatrix (p+1,p+1-freq, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_matrix_pointer += (p+1)*(p+1-freq);

    cublasSetMatrix (p+1-freq, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, p+1-freq, devPtrB, p+1-freq);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_input_pointer +=num_legendre_inputs*(p+1-freq);

    // actual computation:
    cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
	p+1-freq, 1.0F, devPtrA, p+1,
	devPtrB, p+1-freq, 0.0F,
	devPtrC, p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

    cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, &legendre_outputs[leg_outputs_pointer], p+1);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    leg_outputs_pointer += (p+1)*num_legendre_inputs;

    cublasFree (devPtrA);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrB);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
    cublasFree (devPtrC);
    assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  }

  // -------- process  last cosine --------
  num_legendre_inputs = num_vesicles;

  // A will be transformation matrix,  B will be set of inputs (tiled horizontally), C will be set of outputs
  cublasAlloc ((p+1), sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (num_legendre_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc ((p+1)*num_legendre_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasSetMatrix (p+1, 1, sizeof(scalar), (void*)&legendre_matrices[leg_matrix_pointer], p+1, (void*)devPtrA, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_matrix_pointer += (p+1);

  int tmp = (p+1)*(p+1)*(p+2);
  tmp = tmp/2;
  assert(leg_matrix_pointer == tmp);

  cublasSetMatrix (1, num_legendre_inputs, sizeof(scalar), legendre_inputs+leg_input_pointer, 1, devPtrB, 1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_input_pointer +=num_legendre_inputs;
  // assert(leg_input_pointer == 2*p*(p+1)*num_vesicles);
  assert(leg_input_pointer == p*(p+2)*num_vesicles);

  // actual computation:
  cublasSgemm ('N', 'N', p+1, num_legendre_inputs,
      1, 1.0F, devPtrA, p+1,
      devPtrB, 1, 0.0F,
      devPtrC, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (p+1, num_legendre_inputs, sizeof(scalar), devPtrC, p+1, legendre_outputs+leg_outputs_pointer, p+1);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  leg_outputs_pointer += (p+1)*num_legendre_inputs;
  assert(leg_outputs_pointer == 2*p*(p+1)*num_vesicles);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


//   for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//     cout<<legendre_outputs[i]<<endl;

 //  ----- rearrange data for DFT -----
  scalar * dft_inputs = &temp_data2[0];
  leg_input_pointer = 0;

  int vesicle_size = 2*p*(p+1);
  int dft_size = 2*p;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
	dft_inputs[vesicle_size*vesicle_num + dft_size*i + freq] = legendre_outputs[leg_input_pointer++];

  assert(leg_input_pointer==2*p*(p+1)*num_vesicles);

//    for (int i=0; i<2*p*(p+1)*num_vesicles; i++)
//      cout<<dft_inputs[i]<<endl;

  // ---- end of rearranging data ----

  // Now compute backward dft transforms (along each fixed latitude)
  // we use COLUMN-wise storage for dft-matrix

  int num_dft_inputs = num_vesicles*(p+1);
  scalar * dft_outputs = outputs;


  // A will be transformation matrix,  B will be set of inputs (tiled horizontallY), C will be set of dft_outputs
  cublasAlloc (dft_size*dft_size, sizeof(scalar), (void**)&devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasAlloc (dft_size*num_dft_inputs, sizeof(scalar), (void**)&devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  scalar * dft_matrix = &dft_d1backw[0];

  cublasSetMatrix (dft_size, dft_size, sizeof(scalar), (void*)dft_matrix, dft_size, (void*)devPtrA, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasSetMatrix (dft_size, num_dft_inputs, sizeof(scalar), dft_inputs, dft_size, devPtrB, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  // actual computation:
  cublasSgemm ('T', 'N', dft_size, num_dft_inputs,
      dft_size, 1.0F, devPtrA, dft_size,
      devPtrB, dft_size, 0.0F,
      devPtrC, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);


  cublasGetMatrix (dft_size, num_dft_inputs, sizeof(scalar), devPtrC, dft_size, dft_outputs, dft_size);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);

  cublasFree (devPtrA);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrB);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
  cublasFree (devPtrC);
  assert(cublasGetError ()==CUBLAS_STATUS_SUCCESS);
};

