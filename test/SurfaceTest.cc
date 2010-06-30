#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Surface.h"
#include "Parameters.h"

using namespace std;

extern const Device<GPU> the_device(0);
typedef float real;

int main(int argc, char ** argv)
{
    {
        typedef typename containers::Scalars<real,GPU,the_device> Sca;
        typedef typename containers::Vectors<real,GPU,the_device> Vec;
        typedef OperatorsMats<real, Device<GPU> > Opts;
        int const p(12);
        int const nVec(2);
        
        ///@todo Parameters have nothing to do with the surface
        Parameters<real> sim_par;

        //IO
        DataIO<real, Device<GPU> > myIO(the_device);
        
        // initializing vesicle positions from text file
        Vec x0(nVec, p);
        int fLen = x0.getStride();
        
        char fname[400];
        sprintf(fname, "%s/precomputed/dumbbell_cart12_single.txt",
            getenv("VES3D_DIR"));
        myIO.ReadData(fname, fLen * x0.getTheDim(), x0.begin());
        
        for(int ii=1;ii<nVec; ++ii)
            x0.getDevice().Memcpy(x0.getSubN(ii), x0.begin(), x0.getTheDim() * 
                fLen * sizeof(real), MemcpyDeviceToDevice);
        
        //Reading operators from file
        bool readFromFile = true;
        Opts mats(the_device, myIO, readFromFile, sim_par);
        
        //Creating objects
        Surface<Sca, Vec> S(x0, mats);
        Sca X(nVec,p);
        Sca Y(nVec,p);
        Sca Z(nVec,p);

        for(int ii=0;ii<nVec;ii++)
        {
            x0.getDevice().Memcpy(X.getSubN(ii), x0.getSubN(ii), 
                fLen * sizeof(real), MemcpyDeviceToDevice);
            
            x0.getDevice().Memcpy(Y.getSubN(ii), x0.getSubN(ii) + fLen, 
                fLen * sizeof(real), MemcpyDeviceToDevice);
            
            x0.getDevice().Memcpy(Z.getSubN(ii), x0.getSubN(ii) + 2*fLen, 
                fLen * sizeof(real), MemcpyDeviceToDevice); 
        }
        
        //Area and volume
        Sca Area(nVec, p, make_pair(1,1));
        S.area(Area);      

        //the gpu integrator should be fixed
//         Sca Vol(nVec, p, make_pair(1,1));
//         S.volume(Vol);
        cout<<" Area = "<<MaxAbs(Area)<<endl;//", Volume = "<<MaxAbs(Vol)<<endl;
        
        //Vec Cntrs(nVec, 0, make_pair(1,1));
        //S.getCenters(Cntrs);
        //cout<<" Centers :\n"<<Cntrs<<endl;
        
        // Checking the grad and div operator
//         Vec grad(nVec,p), lap(nVec,p);
        
//         S.grad(X,grad);
//         S.div(grad,X);
        
//         S.grad(Y,grad);
//         S.div(grad,Y);
            
//         S.grad(Z,grad);
//         S.div(grad,Z);
    
//         for(int ii=0;ii<nVec;ii++)
//         {
//             lap.getDevice().Memcpy(lap.getSubN(ii), X.getSubN(ii), fLen * 
//                 sizeof(real), MemcpyDeviceToDevice);
            
//             lap.getDevice().Memcpy(lap.getSubN(ii) + fLen, Y.getSubN(ii), fLen * 
//                 sizeof(real), MemcpyDeviceToDevice);
            
//             lap.getDevice().Memcpy(lap.getSubN(ii) + 2*fLen, Z.getSubN(ii), fLen * 
//                 sizeof(real), MemcpyDeviceToDevice);
//         }

//         Sca hh(nVec,p);
//         GeometricDot(lap,S.getNormal(),hh);
//         axpy((real) -.5, hh, S.getMeanCurv(),hh);
        
//         cout<<" The error in the surface grad (For the "
//             <<"\n dumbbell .13120 expected - 2/3 filtering) = "
//             <<fixed<<setprecision(5)<<MaxAbs(hh)<<endl;
        
//         Sca div_n(nVec,p);
//         S.div(S.getNormal(), div_n);
//         axpy((real) .5, div_n, S.getMeanCurv(),div_n);
        
//         cout<<" The error in the surface divergence (For the "
//             <<"\n dumbbell .02964 expected - 2/3 filtering) = "
//             <<fixed<<setprecision(5)<<MaxAbs(div_n)<<endl;

//         S.linearizedMeanCurv(S.getPosition(), hh);
//         axpy((real) -1, hh, S.getMeanCurv(), hh);
//         cout<<" Linear curvature operator: "<<MaxAbs(hh)<<endl;
        
//         grad.getDevice().Memcpy(grad.begin(), S.getNormal().begin(), 
//             S.getNormal().size() * sizeof(real), MemcpyDeviceToDevice);
//         S.mapToTangentSpace(grad);
//         cout<<" Map to tangent space: "<<MaxAbs(grad)<<endl;
    }
    
    PROFILEREPORT(SortTime);
    return 0;
}
