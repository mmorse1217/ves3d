#include "Scalars.h"
#include "HelperFuns.h"

#define T double

using namespace std;

template<typename Container>
void Ax(Container &x, Container &ax)
{
    size_t n = x.Size();
    ax[0] = x[0];
    for(int ii(1); ii<n;++ii)
        ax[ii] = ii* ii * x[ii];
}

int main(int argc, char **argv)
{
    Device<CPU> dev;
    int p = 12;
    int num_funs(1);

    Scalars<T,CPU> x(&dev, p, num_funs), b(&dev, p, num_funs), b2(&dev, p, num_funs);
    int max_iter = 1000;
    T tol = 1e-10;

    for(int ii(0);ii<x.Size();++ii)
        x[ii] = drand48();
    Ax(x,b);

    for(int ii(0);ii<x.Size();++ii)
        x[ii] = 0;

    cout<<BiCGStab(&Ax, x, b, max_iter, tol)<<endl;
    Ax(x,b2);
    axpy((T) -1.0, b, b2, b2);
    
    cout<<"Residual: "<<tol<<endl;
    cout<<"Iter    : "<<max_iter<<endl;
    return 0;
}
