#define get_seconds()   (gettimeofday(&tp, &tzp), \
    (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)

#include <vector>
#include <sys/time.h>
using namespace std;

#define scalar float

class blas_sht
{
  public:
    void forward(scalar * inputs, scalar * outputs);
// "inputs" is array of physical coordinates, vesicle after vesicle;
// coordinates of each vesicle are stored as ( 2p -by- p+1 ) matrix in
// columnwise order. That is, first DFT will be performed on columns
// of each matrix.  "outputs" is array of coefficients in the basis of
// spherical harmonics. Layout: first go Legendre coefs of first
// frequency (constant) for all vesicles, then Legendre coefs of the
// "cos x" term for all vesicles, then "sin x", then "cos 2x" and so
// forth. For each frequency, set of coefficients for the first
// vesicle is followed by set of coefficients for second vesicle, and
// so forth. Note that first frequency has p+1 coefficients, 2nd and
// 3rd have p coefficients, last frequency (sawtooth) has only one
// coefficient.  NOTES: The transform is real (not complex). The
// normalization constant used in DFTs is such that DFT output should
// be multiplied by sines/cosines without further normalization.



    void backward(scalar * inputs, scalar * outputs);
// "inputs" is array of spherical harmonics coefficients.  Layout:
// first go Legendre coefs of first frequency (constant) for all
// vesicles, then Legendre coefs of the "cos \phi" term for all
// vesicles, then "sin \phi", then "cos 2\phi" and so forth. For each
// frequency, set of coefficients for the first vesicle is followed by
// set of coefficients for second vesicle, and so forth. Note that
// first frequency has p+1 coefficients, 2nd and 3rd have p
// coefficients, last frequency (sawtooth) has only one coefficient.
// "outputs" is array of physical coordinates, vesicle after vesicle;
// coordinates of each vesicle are stored as ( 2p -by- p+1 ) matrix in
// columnwise order (different columns correspond to different
// longtitudes).



//     void backward_du(scalar * inputs, scalar * outputs);
    // same for F_u

//     void backward_dv(scalar * inputs, scalar * outputs);
    // F_v
    
//     void backward_d2u(scalar * inputs, scalar * outputs);
    // same for F_uu

//     void backward_d2v(scalar * inputs, scalar * outputs);
    // same for F_vv
    
//     void backward_duv(scalar * inputs, scalar * outputs);
    // same for F_uv

    // constructor allocates temp. memory on host, loads Leg. transforms, calls cublasInit(),  etc.
    blas_sht(int p, int num_vesicles, char * legTransFname, char* legTransInvFname, char* d1legTransFname, char* d2legTransFname);

  private:
    int p;
    int num_vesicles;
    struct timeval  tp;
    struct timezone tzp;
    vector<scalar> legTrans; 
    vector<scalar> legTransInv; 
    vector<scalar> d1legTrans; 
    vector<scalar> d2legTrans; 
    vector<scalar> dft_forward; 
    vector<scalar> dft_backw; 
    vector<scalar> dft_d1backw; 
    vector<scalar> dft_d2backw; 

    vector<scalar>  temp_data;
    vector<scalar>  temp_data2;
};

