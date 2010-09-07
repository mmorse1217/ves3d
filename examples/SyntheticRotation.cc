#include "Logger.h"
#include "Device.h"

#define real_t float
#define DT GPU

void fillRand(real_t* x, int size, const Device<DT> &device)
{
  real_t* y = ( DT == CPU ) ? x : (real_t*) malloc( size * sizeof(real_t));
    
  for(int idx=0;idx<size; ++idx)   
    y[idx] = static_cast<real_t>(drand48());
    
  if ( DT != CPU )
    {
      device.Memcpy(x, y, size * sizeof(real_t), MemcpyHostToDevice);
      free(y);
    }
}

/** 
 * Circular shift of a collection of arrays. The length of x_in is
 * n_sub * sub_length.
 * 
 * @param device A reference to the device on which the arrays reside
 * @param x_in The collection of arrays, stored one after each other
 * @param n_sub Total number of arrays in the memory
 * @param sub_length The length of each array
 * @param shift The shift (of each array)
 * @param x_out The output
 */
void CircShift(const Device<DT> &device, const real_t *x_in, int n_sub, 
	       int sub_length, int shift, real_t *x_out)
{
  PROFILESTART();

  shift = shift % sub_length;
  shift += (shift < 0) ?  sub_length : 0;

  int in_idx, out_idx;
  for (int ii = 0; ii < n_sub; ii++) {
    out_idx = ii * sub_length;
    in_idx = out_idx + sub_length - shift;
    device.Memcpy(x_out + out_idx, x_in + in_idx, sizeof(real_t) * shift, 
		  MemcpyDeviceToDevice);
        
    in_idx = out_idx;
    out_idx += shift;

    device.Memcpy(x_out + out_idx, x_in + in_idx, sizeof(real_t) *
        (sub_length - shift), MemcpyDeviceToDevice);
  }
  PROFILEEND("",0);
}

///The same as above but using transpose
void CircShiftTrans(const Device<DT> &device, const real_t *x_in, int n_sub, 
		    int sub_length, int shift, real_t *x_out)
{
  PROFILESTART();

  shift = shift % sub_length;
  shift += (shift < 0) ?  sub_length : 0;

  device.Transpose(x_in, n_sub, sub_length, x_out);
  real_t* wrk = (real_t*) device.Malloc(n_sub * sub_length * sizeof(real_t));

  int offset = n_sub * (sub_length - shift);
  int allshift = shift * n_sub;

  device.Memcpy(wrk           , x_out + offset, sizeof(real_t) * allshift, 
      MemcpyDeviceToDevice);
  device.Memcpy(wrk + allshift, x_out         , sizeof(real_t) * offset  , 
      MemcpyDeviceToDevice);
    
  device.Transpose(wrk, sub_length, n_sub, x_out);
  PROFILEEND("",0);
}

int main(int , char** )
{
  Device<DT> device(0);     //the device with id = 0, DT is either GPU or CPU

  for ( int p(6); p<24; ++p )   //problem size parameter, 
  {
      cout<<"  Size :"<<p<<endl;
      //Size of the matrices -- blas convention
      int m(2 * p * (p + 1));  
      int n(1024);             //varies between 800~4000
      int k(m);
    
      //Derived sizes
      int a_size(m * k);
      int b_size(k * n);
      int c_size(k * n);
      int all_mats_size((p + 1) * a_size);
 
      //Allocating memory
      real_t *all_mats = (real_t*) device.Malloc(all_mats_size * sizeof(real_t));
      real_t *A        = (real_t*) device.Malloc(a_size        * sizeof(real_t));
      real_t *B        = (real_t*) device.Malloc(b_size        * sizeof(real_t));
      real_t *C        = (real_t*) device.Malloc(c_size        * sizeof(real_t));
    
      //Filling with random numbers
      fillRand(all_mats, all_mats_size, device);
      fillRand(B       , b_size       , device);
    
      //Example of the operation 
      real_t alpha(1), beta(0);
      PROFILECLEAR();
      PROFILESTART();
      for ( int ii=0; ii<p+1; ++ii )
        for ( int jj=0; jj< 2 * p; ++jj)
	  {
            //Permuting the matrix A_i, for the jth entry. Another
            //option is using the transpose based function
            //(CircShiftTrans), but it is msize of the arrays corresponding to the problem sizeuch slower
            CircShift(device, all_mats + ii * a_size, p + 1, 2 * p * k, jj * k, A); 
            device.gemm("N", "N", &m, &n, &k, &alpha, A, &k, B, &k, &beta, C, &k);
	  }
    
      PROFILEEND("",0);
      PROFILEREPORT(SortTime);
    
      //Freeing memory
      device.Free(all_mats);
      device.Free(A);
      device.Free(B);
      device.Free(C);
    }
}
