#include "Scalars.h"
#include "HelperFuns.h"
#include "BiCGStab.h"

#define T double

using namespace std;

#ifndef Doxygen_skip

template<typename Container>
class MatVec
{
 public:
    void operator()(Container &x, Container &ax) const
    {
        size_t n = x.size();
        ax[0] = x[0];
        for(int ii(1); ii<n;++ii)
            ax[ii] = ii* ii * x[ii];
    }
};
#endif //Doxygen_skip

extern const Device<CPU> the_cpu_dev(0);

int main(int argc, char **argv)
{
    typedef containers::Scalars<T,CPU, the_cpu_dev> Sca;
 
    int p = 12;
    int num_funs(1);
    Sca x(num_funs, p), b(num_funs,p), b2(num_funs,p);
    int max_iter = 1000;
    T tol = 1e-10;
    MatVec<Sca> Ax;

    for(int ii(0);ii<x.size();++ii)
        x[ii] = drand48();
    Ax(x,b);

    for(int ii(0);ii<x.size();++ii)
        x[ii] = 0;

    BiCGStab<Sca, MatVec<Sca> > Solver;
    cout<<Solver(Ax, x, b, max_iter, tol)<<endl;

    Ax(x,b2);
    axpy((T) -1.0, b, b2, b2);
    
    cout<<"Residual: "<<tol<<endl;
    cout<<"Iter    : "<<max_iter<<endl;
    return 0;
}

