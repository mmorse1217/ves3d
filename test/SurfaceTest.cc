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

template<typename ST, typename MT>
void test_streaming(const ST &S, const MT &mats){
    std::stringstream s1,s2;
    S.pack(s1, Streamable::ASCII);
    ST S2(S.getShOrder(), mats, NULL,  S.diffFilterFreq(), S.reparamFilterFreq());
    S2.unpack(s1, Streamable::ASCII);
    S2.pack(s2, Streamable::ASCII);
    ASSERT(s1.str()==s2.str(), "bad streaming");
}

template<typename ST>
void test_geo_props(const ST &S){

    typedef typename ST::Sca_t Sca_t;
    typedef typename ST::Vec_t Vec_t;
    typedef typename ST::device_type Dev_t;

    const Vec_t &x0(S.getPosition());
    int p(S.getShOrder());
    int nVec(x0.getNumSubs());
    int fLen(x0.getStride());

    COUT("Extracting X, Y, and Z coordindate functions");
    Sca_t X(nVec,p);
    Sca_t Y(nVec,p);
    Sca_t Z(nVec,p);

    for(int ii=0;ii<nVec;ii++)
    {
        x0.getDevice().Memcpy(X.getSubN_begin(ii),
            x0.getSubN_begin(ii),
            fLen * sizeof(real),
            Dev_t::MemcpyDeviceToDevice);

        x0.getDevice().Memcpy(Y.getSubN_begin(ii),
            x0.getSubN_begin(ii) + fLen,
            fLen * sizeof(real),
            Dev_t::MemcpyDeviceToDevice);

        x0.getDevice().Memcpy(Z.getSubN_begin(ii),
            x0.getSubN_begin(ii) + 2*fLen,
            fLen * sizeof(real),
            Dev_t::MemcpyDeviceToDevice);
    }

    //Area and volume
    real err;
    COUT("Computing area");
    Sca_t Area(nVec, p, std::make_pair(1,1));
    S.area(Area);
    real area(MaxAbs(Area));
    COUT("Area = "<<area);
    ASSERT( fabs(area/16.2179377312307-1)<1e-8,"Expected area for dumbell");

    COUT("Computing volume");
    Sca_t Vol(nVec, p, std::make_pair(1,1));
    S.volume(Vol);
    real vol(MaxAbs(Vol));
    COUT("Volume = "<<vol);
    ASSERT( fabs(vol/5.24886489292959-1)<1e-8,"Expected volume for dumbell");

    COUT("Computing centers");
    Vec_t Cntrs(nVec, 0, std::make_pair(1,1));
    S.getCenters(Cntrs);
    real cntr(MaxAbs(Cntrs));
    ASSERT( fabs(cntr)<1e-8,"Expected center");

    COUT("Computing mean curvature");
    Sca_t H(nVec,p);
    axpy((real) 0, H, S.getMeanCurv(),H);
    err = MaxAbs(H);
    ASSERT(fabs(err/1.376627062-1)<8e-8,"Expected curvature ");

    // Checking the grad and div operator
    Vec_t grad(nVec,p), lap(nVec,p);

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
            Dev_t::MemcpyDeviceToDevice);

        lap.getDevice().Memcpy(lap.getSubN_begin(ii) + fLen,
            Y.getSubN_begin(ii), fLen *
            sizeof(real),
            Dev_t::MemcpyDeviceToDevice);

        lap.getDevice().Memcpy(lap.getSubN_begin(ii) + 2*fLen,
            Z.getSubN_begin(ii), fLen *
            sizeof(real),
            Dev_t::MemcpyDeviceToDevice);
    }

    COUT("Comparing surface Laplacian with curvature");
    Sca_t hh(nVec,p);
    GeometricDot(lap,S.getNormal(),hh);
    err=fabs(MaxAbs(hh)/2.428980517523748-1);
    ASSERT(err<1e-7,"dot(Lap,N)="<<err);

    axpy((real) -.5, hh, S.getMeanCurv(),hh);
    err=fabs(MaxAbs(hh)/0.348685011112687-1);
    ASSERT(err<1e-6,"H-.5*dot(Lap,N)"<<err);

    COUT("Computing Div(N)");
    Sca_t div_n(nVec,p);
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
        Dev_t::MemcpyDeviceToDevice);
    S.mapToTangentSpace(grad);
    ASSERT(MaxAbs(grad)<1e-14,"Normal map to tangent space");
}

template<typename ST, typename MT>
void test_resample(const ST &S, const MT &mats){
    INFO("testing up-sampling");
    ST *S1(NULL), *S2(NULL), *S3(NULL);;
    int p1(mats.p_up_), p2(mats.p_up_+1), p3(mats.p_);
    CHK(S.resample(p1, &S1));
    ASSERT(S1->getShOrder()==p1, "wrong sh order");

    ASSERT(S.resample(p2, &S2) == ErrorEvent::NotImplementedError, "expected error");
    CHK(S1->resample(p3, &S3));

    // typename ST::Vec_t err;
    // err.replicate(S.getPosition());

    // axpy(-1.0, S.getPosition(), S3->getPosition(), err);
    // int stride(err.getStride());
    // for (int i = 0; i < err.getNumSubFuncs(); ++i){
    //     for (int j = 0; j < stride; ++j)
    // 	{
    // 	    double e(fabs(err.begin()[i*stride+j]));
    // 	    std::cout<< ((e < 1e-14) ? 0 : e) <<", ";
    // 	}
    // 	std::cout<<"\n----------"<<std::endl;
    // }
    delete S1;
    delete S2;
    delete S3;

    INFO("testing down-sampling");
}

template<typename Sca, typename Vec, typename Device>
void testSurfaceClass(const Device &dev)
{
    typedef typename Sca::array_type Arr_t;
    typedef OperatorsMats<Arr_t> Mats_t;
    int const nVec(2);

    //IO
    DataIO myIO;

    //@todo Parmeters have nothing to do with the surface
    Parameters<real> sim_par;
    sim_par.sh_order = 6;
    sim_par.upsample_freq = 12;
    ASSERT(sim_par.sh_order==6,"Test only works for p=6");

    // initializing vesicle positions from text file
    Vec x0(nVec, sim_par.sh_order);
    int fLen = x0.getStride();

    char fname[400];
    COUT("Loading initial shape");
    sprintf(fname, "precomputed/dumbbell_%d_double.txt",sim_par.sh_order);
    myIO.ReadData(FullPath(fname), x0, DataIO::ASCII, 0, fLen * DIM);

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
    Surface<Sca, Vec> S(x0.getShOrder(),mats, &x0);

    test_geo_props(S);
    test_streaming(S, mats);
    test_resample(S, mats);
}

int main(int argc, char ** argv)
{
    COUT("Surface test:\n=============");
    COUT("CPU device:\n------------");

    typedef Scalars<real, DevCPU, the_cpu_device> ScaCPU_t;
    typedef Vectors<real, DevCPU, the_cpu_device> VecCPU_t;

    testSurfaceClass<ScaCPU_t, VecCPU_t, DevCPU>(the_cpu_device);
    PROFILEREPORT(SortTime);
    COUT(emph<<" *** Surface class with CPU device passed ***"<<emph);

#ifdef GPU_ACTIVE
    PROFILECLEAR();
    COUT("GPU device:\n------------");
    typedef Scalars<real, DevGPU, the_gpu_device> ScaGPU_t;
    typedef Vectors<real, DevGPU, the_gpu_device> VecGPU_t;

    testSurfaceClass<ScaGPU_t, VecGPU_t, DevGPU>(the_gpu_device);
    PROFILEREPORT(SortTime);

    COUT(emph<<" *** Surface class with GPU device passed ***"<<emph);
#endif //GPU_ACTIVE

    return 0;
}
