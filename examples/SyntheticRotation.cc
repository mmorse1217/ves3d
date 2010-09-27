#include "Logger.h"
#include "Device.h"

#define real_t float
#define DT CPU

///Fills the given array with random numbers
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

/** 
 * Permutes the given matrix Ai to the desired form Aij
 * 
 * @param device The device on which the data resides
 * @param Ai The source matrix
 * @param p The size argument
 * @param j The index of the desired permutation
 * @param Aij The output
 */
void Permute(const Device<DT> &device, const real_t *Ai, int p, int j, real_t *Aij)
{
    int dim( 2 * p *(p+1));
    int stride(2*p);
    //The matrices is assumed to be stored column-wise
    //we are exchanging columns (modulus the strides)

    //Going through columns
    for(int jj=0;jj<dim;++jj)
    {
        //Finding the column number in the original matrix
        int original_jj = jj%stride - j;
        original_jj += (original_jj < 0) ? stride : 0;
        original_jj += jj/stride * stride;
             
        //going through rows
        device.Memcpy(Aij + jj * dim, Ai + original_jj * dim , sizeof(real_t) *
            dim, MemcpyDeviceToDevice);
    }
}

int main(int , char** )
{
    Device<DT> device(0);     //the device with id = 0, DT is either GPU or CPU

    int p_min(6), p_max(20);
    for ( int p=p_min; p<p_max; ++p )   //problem size parameter, 
    {
        cout<<"  Size :"<<p<<endl;
        
        //Size of the matrices -- blas convention
        int m(2 * p * (p + 1));  
        int n(1);             //varies between 800~4000
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
                //Permuting the matrix A_i, for the jth entry. 
                
                //CircShift and Permute perform the same task, but
                //Permute is a little more explicit in its operations
                //and easier to understand
                
                //CircShift(device, all_mats + ii * a_size, p + 1, 2 * p * k, jj * k, A);
                Permute(device, all_mats + ii * a_size, p, jj, A);
                
                device.gemm("N", "N", &m, &n, &k, &alpha, A, &k, B, &k, &beta, C, &k);
                
                //Here we would do some processing on C
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
