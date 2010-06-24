#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "Surface.h"
#include "Parameters.h"

using namespace std;

extern const Device<CPU> the_cpu_dev(0);

typedef double real;

int main(int argc, char ** argv)
{
    {
        typedef containers::Scalars<real,CPU,the_cpu_dev>::value_type T;
        typedef typename containers::Scalars<real,CPU,the_cpu_dev> Sca;
        typedef typename containers::Vectors<real,CPU,the_cpu_dev> Vec;

        int const p(12);
        int const nVec(2);
        
        //IO
        DataIO<T,CPU> myIO(the_cpu_dev);
        
        // initializing vesicle positions from text file
        Vec x0(nVec, p);
        int fLen = x0.getStride();
        
        char fname[400];
        sprintf(fname, "%s/precomputed/dumbbell_cart12_single.txt",
            getenv("VES3D_DIR"));
        myIO.ReadData(fname, fLen * x0.getTheDim(), x0.begin());
        
        for(int ii=1;ii<nVec; ++ii)
            x0.getDevice().Memcpy(x0.getSubN(ii), x0.begin(), x0.getTheDim() * 
                fLen * sizeof(T), MemcpyDeviceToDevice);
        
        //Reading operators from file
        bool readFromFile = true;
        OperatorsMats<real> mats(myIO, readFromFile);
        
        //Creating objects
        Surface<Sca, Vec> S(x0, mats);
        Sca X(nVec,p);
        Sca Y(nVec,p);
        Sca Z(nVec,p);

        for(int ii=0;ii<nVec;ii++)
        {
            x0.getDevice().Memcpy(X.getSubN(ii), x0.getSubN(ii), 
                fLen * sizeof(T), MemcpyDeviceToDevice);
            
            x0.getDevice().Memcpy(Y.getSubN(ii), x0.getSubN(ii) + fLen, 
                fLen * sizeof(T), MemcpyDeviceToDevice);
            
            x0.getDevice().Memcpy(Z.getSubN(ii), x0.getSubN(ii) + 2*fLen, 
                fLen * sizeof(T), MemcpyDeviceToDevice);
        }
        
        //Area and volume
        Sca Area(nVec, p, make_pair(1,1));
        S.area(Area);
        cout<<" Area : "<<Area[0]<<endl;
        
        Sca Vol(nVec, p, make_pair(1,1));
        S.volume(Vol);
        cout<<" Volume : "<<Vol[0]<<endl;
        
        Vec Cntrs(nVec, 0, make_pair(1,1));
        S.getCenters(Cntrs);
        //cout<<" Centers :\n"<<Cntrs<<endl;
        
        // Checking the grad and div operator
        Vec grad(nVec,p), lap(nVec,p);
        
        S.grad(X,grad);
        S.div(grad,X);
        
        S.grad(Y,grad);
        S.div(grad,Y);
            
        S.grad(Z,grad);
        S.div(grad,Z);
    
        for(int ii=0;ii<nVec;ii++)
        {
            lap.getDevice().Memcpy(lap.getSubN(ii), X.getSubN(ii), fLen * 
                sizeof(T), MemcpyDeviceToDevice);
            
            lap.getDevice().Memcpy(lap.getSubN(ii) + fLen, Y.getSubN(ii), fLen * 
                sizeof(T), MemcpyDeviceToDevice);
            
            lap.getDevice().Memcpy(lap.getSubN(ii) + 2*fLen, Z.getSubN(ii), fLen * 
                sizeof(T), MemcpyDeviceToDevice);
        }

        Sca hh(nVec,p);
        GeometricDot(lap,S.getNormal(),hh);
        axpy((T) -.5, hh, S.getMeanCurv(),hh);
        
        cout<<"\n The error in the surface grad (For the dumbbell"
            <<" .13120 expected - 2/3 filtering)= "
            <<hh.getDevice().MaxAbs(hh.begin(), hh.size())<<endl;
 
        Sca div_n(nVec,p);
        S.div(S.getNormal(), div_n);
        axpy((T) .5, div_n, S.getMeanCurv(),div_n);
        
        cout<<"\n The error in the surface divergence (For the "
            <<"dumbbell .02964 expected - 2/3 filtering)= "
            <<div_n.getDevice().MaxAbs(div_n.begin(), div_n.size())<<endl;

        S.linearizedMeanCurvature(S.getPosition(), hh);
        axpy((T) -1, hh, S.getMeanCurv(), hh);
    }

    //PROFILEREPORT(SortTime);
    return 0;
}
