#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"

using namespace std;

extern const Device<CPU> the_cpu_dev(0);

int main(int argc, char ** argv)
{
    {
        typedef containers::Scalars<float,CPU,the_cpu_dev>::value_type T;
        typedef typename containers::Scalars<float,CPU,the_cpu_dev> Sca;
        typedef typename containers::Vectors<float,CPU,the_cpu_dev> Vec;

        int const p(12);
        int const nVec(1);
        
        //IO
        DataIO<T,CPU> myIO(the_cpu_dev);
        
        // initializing vesicle positions from text file
        Vec x0(nVec, p);
                       
        char fname[400];
        //sprintf(fname,"%s/precomputed/sphere_cart12_single.txt",getenv("VES3D_DIR"));
        //sprintf(fname,"%s/precomputed/two_ellipse_12",getenv("VES3D_DIR"));
        sprintf(fname,"%s/precomputed/dumbbell_cart12_single.txt",getenv("VES3D_DIR"));
        myIO.ReadData(fname,x0.size(),x0.begin());

        //Creating objects
        Surface<Sca, Vec> S(x0);
        Sca X(nVec,p);
        Sca Y(nVec,p);
        Sca Z(nVec,p);

        size_t fLen = x0.getStride();
        for(int ii=0;ii<nVec;ii++)
            for(int idx=0;idx<fLen;idx++)
            {
                *(X.begin() +  ii*fLen + idx) = *(x0.begin() + idx);
                *(Y.begin() +  ii*fLen + idx) = *(x0.begin() + idx + fLen);
                *(Z.begin() +  ii*fLen + idx) = *(x0.begin() + idx + fLen + fLen);
            }
        
        //Area and volume
        Sca Area(nVec, p, make_pair(1,1));
        S.area(Area);
        //cout<<" Area :\n"<<Area<<endl;
        
        Sca Vol(nVec, p, make_pair(1,1));
        S.volume(Vol);
        //cout<<" Volume :\n"<<Vol<<endl;
        
        Vec Cntrs(nVec, 0, make_pair(1,1));
        S.getCenters(Cntrs);
        //cout<<" Centers :\n"<<Cntrs<<endl;
        
        T err = 0, err2;
        // Checking the grad and div operator
        Vec grad(nVec,p), lap(nVec,p);

        S.grad(X,grad);
        S.div(grad,X);
    
        S.grad(Y,grad);
        S.div(grad,Y);
    
        S.grad(Z,grad);
        S.div(grad,Z);
    
        for(int ii=0;ii<nVec;ii++)
            for(int idx=0;idx<fLen;idx++)
            {
                *(lap.begin() +  3*ii*fLen + idx              ) = *(X.begin() + ii*fLen + idx);
                *(lap.begin() +  3*ii*fLen + idx + fLen       ) = *(Y.begin() + ii*fLen + idx);
                *(lap.begin() +  3*ii*fLen + idx + fLen + fLen) = *(Z.begin() + ii*fLen + idx);
            }

        Sca hh(nVec,p);
        GeometricDot(lap,S.getNormal(),hh);
        axpy((T) -.5, hh, S.getMeanCurv(),hh);
        
        err = 0;
        for(int ii=0;ii<hh.size(); ii++)
        {
            err2 = fabs(hh.begin()[ii]);
            err = (err>err2) ? err : err2;
        }
        cout<<"\n The error in the surface grad (For the dumbbell .13703 expected)= "<<err<<endl;
 
        Sca div_n(nVec,p);
        S.div(S.getNormal(), div_n);
        axpy((T) .5, div_n, S.getMeanCurv(),div_n);
        
        err = 0;
        for(int ii=0;ii<div_n.size(); ii++)
        {
            err2 = fabs(*(div_n.begin() + ii));
            err = (err>err2) ? err : err2;
        }

        cout<<"\n The error in the surface divergence (For the dumbbell .02964 expected)= "<<err<<endl;
    }
    PROFILEREPORT(SortTime);
    return 0;
}
