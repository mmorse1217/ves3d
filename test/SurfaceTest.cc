#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Surface.h"
#include "Parameters.h"

typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);

#ifdef GPU_ACTIVE
typedef Device<GPU> DevGPU;
extern const DevGPU the_gpu_device(0);
#endif //GPU_ACTIVE

typedef double real;

#ifndef Doxygen_skip

template<typename Sca, typename Vec, typename Device>
void testSurface(const Device &dev)
{
    typedef typename Sca::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    int const nVec(2);

    //IO
    DataIO myIO;

    //@todo Parmeters have nothing to do with the surface
    Parameters<real> sim_par;
    sim_par.sh_order = 6;
    sim_par.rep_up_freq = 12;
    ASSERT(sim_par.sh_order==6,"Test only works for p=6");

    // initializing vesicle positions from text file
    Vec x0(nVec, sim_par.sh_order);
    int fLen = x0.getStride();

    char fname[400];
    COUT("Loading initial shape");
    sprintf(fname, "%s/precomputed/dumbbell_%d_double.txt",VES3D_PATH,sim_par.sh_order);
    myIO.ReadData(fname, x0, DataIO::ASCII, 0, fLen * DIM);

    COUT("Populating x0 (nves="<<nVec<<")");
    for(int ii=1;ii<nVec; ++ii)
        x0.getDevice().Memcpy(x0.getSubN_begin(ii),
            x0.begin(),
            x0.getTheDim() * fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);

    //Reading operators from file
    COUT("Loading matrices");
    bool readFromFile = true;
    Mats_t mats(readFromFile, sim_par);

    //Creating objects
    COUT("Creating the surface object");
    Surface<Sca, Vec> S(x0, mats);

    COUT("Extracting X, Y, and Z coordindate functions");
    Sca X(nVec,sim_par.sh_order);
    Sca Y(nVec,sim_par.sh_order);
    Sca Z(nVec,sim_par.sh_order);

    for(int ii=0;ii<nVec;ii++)
    {
        x0.getDevice().Memcpy(X.getSubN_begin(ii),
            x0.getSubN_begin(ii),
            fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);

        x0.getDevice().Memcpy(Y.getSubN_begin(ii),
            x0.getSubN_begin(ii) + fLen,
            fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);

        x0.getDevice().Memcpy(Z.getSubN_begin(ii),
            x0.getSubN_begin(ii) + 2*fLen,
            fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);
     }

    //Area and volume
    real err;
    COUT("Computing area");
    Sca Area(nVec, sim_par.sh_order, std::make_pair(1,1));
    S.area(Area);
    real area(MaxAbs(Area));
    COUT("Area = "<<area);
    ASSERT( fabs(area/16.21793733-1)<1e-8,"Expected area for dumbell");

    COUT("Computing volume");
    Sca Vol(nVec, sim_par.sh_order, std::make_pair(1,1));
    S.volume(Vol);
    real vol(MaxAbs(Vol));
    COUT("Volume = "<<vol);
    ASSERT( fabs(vol/5.24886478864-1)<1e-8,"Expected area for dumbell");

    COUT("Computing centers");
    Vec Cntrs(nVec, 0, std::make_pair(1,1));
    S.getCenters(Cntrs);
    real cntr(MaxAbs(Cntrs));
    ASSERT( fabs(cntr)<1e-8,"Expected center");

    COUT("Computing mean curvature");
    Sca H(nVec,sim_par.sh_order);
    axpy((real) 0, H, S.getMeanCurv(),H);
    err = MaxAbs(H);
    ASSERT(fabs(err/1.376627062-1)<8e-8,"Expected curvature ");

    // Checking the grad and div operator
    Vec grad(nVec,sim_par.sh_order), lap(nVec,sim_par.sh_order);

    COUT("Computing surface Laplacian");
    S.grad(X,grad);
    err=fabs(MaxAbs(grad)/1.01100423438481-1);
    ASSERT(err<1e-7,"grad X error="<<err);
    S.div(grad,X);
    err=fabs(MaxAbs(X)/2.239995450856133-1);
    ASSERT(err<1e-7,"Laplacian X error="<<err);

    S.grad(Y,grad);
    err=fabs(MaxAbs(grad)/1.011004234384816-1);
    ASSERT(err<1e-7,"grad Y error="<<err);
    S.div(grad,Y);
    err=fabs(MaxAbs(Y)/2.239995450856133-1);
    ASSERT(err<1e-7,"Laplacian Y error="<<err);

    S.grad(Z,grad);
    err=fabs(MaxAbs(grad)/1.004308085588217-1);
    ASSERT(err<1e-7,"grad Z error="<<err);
    S.div(grad,Z);
    err=fabs(MaxAbs(Z)/1.777133873450119-1);
    ASSERT(err<1e-7,"Laplacian Z error="<<err);

    for(int ii=0;ii<nVec;ii++)
    {
        lap.getDevice().Memcpy(lap.getSubN_begin(ii),
            X.getSubN_begin(ii),
            fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);

        lap.getDevice().Memcpy(lap.getSubN_begin(ii) + fLen,
            Y.getSubN_begin(ii), fLen *
            sizeof(real),
            Device::MemcpyDeviceToDevice);

        lap.getDevice().Memcpy(lap.getSubN_begin(ii) + 2*fLen,
            Z.getSubN_begin(ii), fLen *
            sizeof(real),
            Device::MemcpyDeviceToDevice);
    }

    COUT("Comparing surface Laplacian with curvature");
    Sca hh(nVec,sim_par.sh_order);
    GeometricDot(lap,S.getNormal(),hh);
    err=fabs(MaxAbs(hh)/2.428980517523748-1);
    ASSERT(err<1e-7,"dot(Lap,N)="<<err);

    axpy((real) -.5, hh, S.getMeanCurv(),hh);
    err=fabs(MaxAbs(hh)/0.348685011112687-1);
    ASSERT(err<1e-6,"H-.5*dot(Lap,N)"<<err);

    COUT("Computing Div(N)");
    Sca div_n(nVec,sim_par.sh_order);
    S.div(S.getNormal(), div_n);
    axpy((real) .5, div_n, S.getMeanCurv(),div_n);
    err = MaxAbs(div_n);
    ASSERT(fabs(err/0.2437253515-1)<1e-6,"Expected error");

    COUT("Checking linearizedMeanCurv");
    S.linearizedMeanCurv(S.getPosition(), hh);
    axpy((real) -1, hh, S.getMeanCurv(), hh);
    ASSERT(fabs(MaxAbs(hh))<1e-14,"linear curvature is the same as H");

    COUT("Checking mapToTangentSpace");
    grad.getDevice().Memcpy(grad.begin(), S.getNormal().begin(),
        S.getNormal().size() * sizeof(real),
        Device::MemcpyDeviceToDevice);
    S.mapToTangentSpace(grad);
    ASSERT(MaxAbs(grad)<1e-14,"Normal map to tangent space");
 }
 #endif //Doxygen_skip

 int main(int argc, char ** argv)
 {
     COUT("Surface test:\n=============");
     COUT("CPU device:\n------------");

     typedef Scalars<real, DevCPU, the_cpu_device> ScaCPU_t;
     typedef Vectors<real, DevCPU, the_cpu_device> VecCPU_t;

    testSurface<ScaCPU_t, VecCPU_t, DevCPU>(the_cpu_device);
    PROFILEREPORT(SortTime);
    COUT(emph<<" *** Surface class with CPU device passed ***"<<emph);

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("GPU device:\n------------");
    typedef Scalars<real, DevGPU, the_gpu_device> ScaGPU_t;
    typedef Vectors<real, DevGPU, the_gpu_device> VecGPU_t;

    testSurface<ScaGPU_t, VecGPU_t, DevGPU>(the_gpu_device);
    PROFILEREPORT(SortTime);

    COUT(emph<<" *** Surface class with GPU device passed ***"<<emph);
#endif //GPU_ACTIVE

    return 0;
}
