#include "Scalars.h"
#include "HelperFuns.h"
#include "BiCGStab.h"
#include <stdlib.h>

typedef float real;

using namespace std;

#ifndef Doxygen_skip

template<typename Container>
class MatVec
{
  public:
    typedef typename Container::value_type T;
    Container diag;

    MatVec(int nf, int p)
    {
        diag.resize(nf, p);
        size_t size = diag.size();

        T* buffer = new T[size];       
        for(int ii=0; ii<size;++ii)
            buffer[ii] = (ii + 1) * (ii + 1);
        
        diag.getDevice().Memcpy(diag.begin(), buffer, 
            size * sizeof(T), MemcpyHostToDevice);
        
        delete[] buffer;
    }
    
    void operator()(Container &x, Container &ax) const
    {
        xy(diag, x, ax);
    }
};
#endif //Doxygen_skip

extern const Device<CPU> the_cpu_dev(0);

#ifdef GPU_ACTIVE
extern const Device<GPU> the_gpu_dev(0);
#endif //GPU_ACTIVE

int main(int argc, char **argv)
{
    COUT("\n ==============================\n"
        <<"  BiCGStab Test:"
        <<"\n ==============================\n");
    sleep(.5);
    
    
    typedef Scalars<real,CPU, the_cpu_dev> ScaCPU;
    
    int p = 12;
    int nfuns(1);
    ScaCPU x_ref(nfuns,p), b_ref(nfuns,p);
    MatVec<ScaCPU> Ax(nfuns, p);
    const int max_iter = 200;
    const real tol = 1e-6;

    for(int ii(0);ii<x_ref.size();++ii)
        *(x_ref.begin() + ii) = drand48();
    Ax(x_ref, b_ref);
    
    char cpu_out[300], gpu_out[300];
    
    {
        ScaCPU x(nfuns, p), b(nfuns,p);
        b.getDevice().Memcpy(b.begin(), b_ref.begin(),
            b.size() * sizeof(real), MemcpyHostToDevice);
        
        BiCGStab<ScaCPU, MatVec<ScaCPU> > Solver;
        int miter(max_iter);
        real tt(tol);
    
        enum BiCGSReturn ret = Solver(Ax, x, b, miter, tt);
       
        Ax(x,b);
        axpy((real) -1.0, b_ref, b, b);
                
        string formatstr ="\n CPU data :\n     Residual: %2.4e\n     ";
        formatstr +="Iter    : %d\n     Error   : %2.4e\n";
        sprintf(cpu_out, formatstr.c_str(), tt, miter, MaxAbs(b));
    }
    
#ifdef GPU_ACTIVE
    {
        typedef containers::Scalars<real,GPU, the_gpu_dev> ScaGPU;
        MatVec<ScaGPU> AxGPU(nfuns, p);
        
        ScaGPU x(nfuns, p), b(nfuns,p);
        b.getDevice().Memcpy(b.begin(), b_ref.begin(),
            b.size() * sizeof(real), MemcpyHostToDevice);

        axpy(0, x, x);

        BiCGStab<ScaGPU, MatVec<ScaGPU> > Solver;
        int miter(max_iter);
        real tt(tol);
        
        enum BiCGSReturn ret = Solver(AxGPU, x, b, miter, tt);
        
        AxGPU(x,b);
        
        ScaCPU b_cpu(nfuns, p);
        b.getDevice().Memcpy(b_cpu.begin(), b.begin(),
            b.size() * sizeof(real), MemcpyDeviceToHost);

        axpy((real) -1.0, b_ref, b_cpu, b_cpu);
        
        string formatstr ="\n GPU data :\n     Residual: %2.4e\n     ";
         formatstr +="Iter    : %d\n     Error   : %2.4e\n";
        sprintf(gpu_out, formatstr.c_str(), tt, miter, MaxAbs(b_cpu));
    }
#endif //GPU_ACTIVE
    
    cout<<cpu_out;
#ifdef GPU_ACTIVE
    cout<<gpu_out;
#endif //GPU_ACTIVE
    cout<<endl;
    sleep(.5);
    
   return 0;
}

