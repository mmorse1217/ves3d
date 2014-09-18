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
extern const Device<GPU> the_gpu_device(0);
#endif //GPU_ACTIVE

typedef double real;

#ifndef Doxygen_skip

template<typename Sca, typename Vec, typename Device>
void testSurface(const Device &dev)
{
    typedef typename Sca::array_type Arr;
    typedef OperatorsMats<Arr> Mats_t;
    int const nVec(2);

    //IO
    DataIO myIO;

    //@todo Parmeters have nothing to do with the surface
    Parameters<real> sim_par;
    sim_par.sh_order = 6;
    sim_par.rep_up_freq = 12;
    COUT(sim_par<<std::endl);

    // initializing vesicle positions from text file
    Vec x0(nVec, sim_par.sh_order);
    int fLen = x0.getStride();

    char fname[400];
    sprintf(fname, "precomputed/biconcave_ra85_%u",sim_par.sh_order);
    myIO.ReadData(fname, x0, 0, fLen * DIM);

    for(int ii=1;ii<nVec; ++ii)
        x0.getDevice().Memcpy(x0.getSubN_begin(ii),
            x0.begin(),
            x0.getTheDim() * fLen * sizeof(real),
            Device::MemcpyDeviceToDevice);

    //Reading operators from file
    bool readFromFile = true;
    Mats_t mats(readFromFile, sim_par);

    //Creating objects
    Surface<Sca, Vec> S(x0, mats);

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
    Sca Area(nVec, sim_par.sh_order, std::make_pair(1,1));
    S.area(Area);

     //the gpu integrator should be fixed
     //         Sca Vol(nVec, p, make_pair(1,1));
     //         S.volume(Vol);
     COUT(" Area = "<<MaxAbs(Area)<<std::endl);//", Volume = "<<MaxAbs(Vol)<<std::endl;

     //Vec Cntrs(nVec, 0, make_pair(1,1));
     //S.getCenters(Cntrs);
     //COUT(" Centers :\n"<<Cntrs<<std::endl);

     // Checking the grad and div operator
     Vec grad(nVec,sim_par.sh_order), lap(nVec,sim_par.sh_order);

     S.grad(X,grad);
     S.div(grad,X);

     S.grad(Y,grad);
     S.div(grad,Y);

     S.grad(Z,grad);
     S.div(grad,Z);

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

     Sca hh(nVec,sim_par.sh_order);
     GeometricDot(lap,S.getNormal(),hh);
     axpy((real) -.5, hh, S.getMeanCurv(),hh);

     COUT(" The error in the surface grad (For the "
         <<"\n dumbbell .13120 expected - 2/3 filtering) = "
         <<std::fixed<<std::setprecision(5)<<MaxAbs(hh)<<std::endl);

     Sca div_n(nVec,sim_par.sh_order);
     S.div(S.getNormal(), div_n);
     axpy((real) .5, div_n, S.getMeanCurv(),div_n);

     COUT(" The error in the surface divergence (For the "
         <<"\n dumbbell .02964 expected - 2/3 filtering) = "
         <<std::fixed<<std::setprecision(5)<<MaxAbs(div_n)<<std::endl);

     S.linearizedMeanCurv(S.getPosition(), hh);
     axpy((real) -1, hh, S.getMeanCurv(), hh);
     COUT(" Linear curvature operator: "<<MaxAbs(hh)<<std::endl);

     grad.getDevice().Memcpy(grad.begin(), S.getNormal().begin(),
         S.getNormal().size() * sizeof(real),
         Device::MemcpyDeviceToDevice);
     S.mapToTangentSpace(grad);
     COUT(" Map to tangent space: "<<MaxAbs(grad)<<std::endl);
 }
 #endif //Doxygen_skip

 int main(int argc, char ** argv)
 {
     COUT("\n\n ================\n  Surface test: \n ================"<<std::endl);
     COUT("\n ------------ \n  CPU device: \n ------------"<<std::endl);

     typedef Scalars<real, DevCPU, the_cpu_device> ScaCPU_t;
     typedef Vectors<real, DevCPU, the_cpu_device> VecCPU_t;

    testSurface<ScaCPU_t, VecCPU_t, DevCPU>(the_cpu_device);

    PROFILEREPORT(SortTime);

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("\n ------------ \n  GPU device: \n ------------"<<std::endl);
    typedef Scalars<real, DevGPU, the_gpu_device> ScaGPU_t;
    typedef Vectors<real, DevGPU, the_gpu_device> VecGPU_t;

    testSurface<ScaGPU_t, VecGPU_t, DevGPU>(the_gpu_device);

    PROFILEREPORT(SortTime);
#endif //GPU_ACTIVE

    return 0;
}
