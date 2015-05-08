#include <iostream>
#include "DataIO.h"
#include "Scalars.h"
#include "Vectors.h"
#include "OperatorsMats.h"
#include "Surface.h"
#include "Parameters.h"
#include "MovePole.h"

#ifdef GPU_ACTIVE
#include "CudaKernels.h"
#endif //GPU_ACTIVE

typedef Device<CPU> DevCPU;
extern const DevCPU the_cpu_device(0);

#ifdef GPU_ACTIVE
typedef Device<GPU> DevGPU;
extern const DevGPU the_gpu_device(0);
#endif //GPU_ACTIVE

typedef double real;

template<typename Sca_t, typename Vec_t, typename Device>
void testStokesDoubleLayer(const Device &dev){

  typedef Surface<Sca_t, Vec_t> Surf_t;
  typedef typename Surf_t::device_type device_type;

  typedef typename Sca_t::array_type Arr_t;
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
  Vec_t x0(nVec, sim_par.sh_order);
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
  if(1){ // Analytical coordinates (Sphere)
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    for(size_t k=0;k<nVec;k++){
      real* x_k=x0.getSubN_begin(k)+0*fLen;
      real* y_k=x0.getSubN_begin(k)+1*fLen;
      real* z_k=x0.getSubN_begin(k)+2*fLen;
      for(size_t i=0;i<imax;i++)
      for(size_t j=0;j<jmax;j++){
        x_k[j+i*jmax]=cos((i+1)*M_PI/(imax+1));
        y_k[j+i*jmax]=sin((i+1)*M_PI/(imax+1))*sin(j*2*M_PI/jmax);
        z_k[j+i*jmax]=sin((i+1)*M_PI/(imax+1))*cos(j*2*M_PI/jmax);
      }
    }
  }

  //Reading operators from file
  COUT("Loading matrices");
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  //Creating objects
  COUT("Creating the surface object");
  Surf_t S(x0.getShOrder(),mats,&x0);

  Sca_t w_sph_, w_sph_inv_, sing_quad_weights_;
  { // Set w_sph_, w_sph_inv_, sing_quad_weights_
    typedef typename Surf_t::value_type value_type;
    int p = S.getPosition().getShOrder();
    int np = S.getPosition().getStride();

    w_sph_.resize(1, p);
    w_sph_inv_.resize(1, p);
    w_sph_.getDevice().Memcpy(w_sph_.begin(), mats.w_sph_,
        np * sizeof(value_type), device_type::MemcpyDeviceToDevice);
    xInv(w_sph_,w_sph_inv_);

    //Singular quadrature weights
    sing_quad_weights_.resize(1,p);
    sing_quad_weights_.getDevice().Memcpy(sing_quad_weights_.begin(),
        mats.sing_quad_weights_, sing_quad_weights_.size() *
        sizeof(value_type),
        device_type::MemcpyDeviceToDevice);
  }

  MovePole<Sca_t,Mats_t> move_pole(mats);

  //===========================================================================

  Vec_t force   (nVec,sim_par.sh_order);
  Vec_t velocity(nVec,sim_par.sh_order);
  Vec_t velocity_exact(nVec,sim_par.sh_order);

  { // Set force.
    real* force_=new real[force.size()];
    for(size_t i=0;i<force.size();i++){
      force_[i]=1.0;
    }
    force.getDevice().Memcpy(force.begin(), force_,
        force.size() * sizeof(real), device_type::MemcpyHostToDevice);
    delete[] force_;
  }
  { // Set velocity_exact.
    real* velocity_exact_=new real[velocity_exact.size()];
    for(size_t i=0;i<velocity_exact.size();i++){
      velocity_exact_[i]=0.5;
    }
    velocity_exact.getDevice().Memcpy(velocity_exact.begin(), velocity_exact_,
        velocity_exact.size() * sizeof(real), device_type::MemcpyHostToDevice);
    delete[] velocity_exact_;
  }

  int imax(S.getPosition().getGridDim().first);
  int jmax(S.getPosition().getGridDim().second);
  int np = S.getPosition().getStride();
  int nv = S.getPosition().getNumSubs();

  Sca_t t1(nVec,sim_par.sh_order);
  Sca_t t2(nVec,sim_par.sh_order);
  Vec_t v1(nVec,sim_par.sh_order);
  Vec_t v2(nVec,sim_par.sh_order);
  Vec_t v3(nVec,sim_par.sh_order);

  ax(w_sph_inv_, S.getAreaElement(), t1);

  int numinputs = 4;
  const Sca_t* inputs[] = {&S.getPosition(), &S.getNormal(), &force, &t1};
  Sca_t*      outputs[] = {&             v1, &           v3, &   v2, &t2};
  move_pole.setOperands(inputs, numinputs, sim_par.singular_stokes);

  for(int ii=0;ii < imax; ++ii)
  for(int jj=0;jj < jmax; ++jj){
    move_pole(ii, jj, outputs);

    ax(w_sph_, t2, t2);
    xv(t2, v2, v2);

    S.getPosition().getDevice().DirectStokesDoubleLayer(v1.begin(), v3.begin(), v2.begin(),
        sing_quad_weights_.begin(), np, nv, S.getPosition().begin(),
        ii * jmax + jj, ii * jmax + jj + 1, velocity.begin());
  }

  { // Compute error.
    Vec_t err_vec(nVec,sim_par.sh_order);
    axpy((real)(-1.0),velocity,velocity_exact,err_vec);
    //operator<<(std::cout,*(Array<real,Device,the_cpu_device>*)&err_vec);
    real err=MaxAbs(err_vec);
    ASSERT(err<0.03,"Velocity error="<<err);
  }
}

int main(int argc, char ** argv){

  COUT("Surface test:\n=============");
  COUT("CPU device:\n------------");

  typedef Scalars<real, DevCPU, the_cpu_device> ScaCPU_t;
  typedef Vectors<real, DevCPU, the_cpu_device> VecCPU_t;

  testStokesDoubleLayer<ScaCPU_t, VecCPU_t, DevCPU>(the_cpu_device);
  PROFILEREPORT(SortTime);
  COUT(emph<<" *** StokesDoubleLayer with CPU device passed ***"<<emph);

#ifdef GPU_ACTIVE
  cudaError_t error=cudaGetLastError();
  if (error != cudaSuccess) fprintf(stderr,"CUDA Error: %s \n", cudaGetErrorString(error));
  assert(error == cudaSuccess);

  PROFILECLEAR();
  COUT("GPU device:\n------------");
  typedef Scalars<real, DevGPU, the_gpu_device> ScaGPU_t;
  typedef Vectors<real, DevGPU, the_gpu_device> VecGPU_t;

  // GPU code does not work because of bug in:
  // SHTrans<Container, Mats>::collectSameOrder
  //testStokesDoubleLayer<ScaGPU_t, VecGPU_t, DevGPU>(the_gpu_device);
  //PROFILEREPORT(SortTime);

  { // Compare CPU and GPU output
    int nVec=2;
    int sh_order=32;
    VecCPU_t trg_cpu(nVec,sh_order); fillRand(trg_cpu);
    VecCPU_t src_cpu(nVec,sh_order); fillRand(src_cpu);
    VecCPU_t den_cpu(nVec,sh_order); fillRand(den_cpu);
    VecCPU_t nor_cpu(nVec,sh_order); fillRand(nor_cpu);
    VecCPU_t  qw_cpu(nVec,sh_order); fillRand( qw_cpu);
    VecCPU_t vel_cpu(nVec,sh_order); //fillRand(vel_cpu);

    VecGPU_t trg_gpu(nVec,sh_order);
    VecGPU_t src_gpu(nVec,sh_order);
    VecGPU_t den_gpu(nVec,sh_order);
    VecGPU_t nor_gpu(nVec,sh_order);
    VecGPU_t  qw_gpu(nVec,sh_order);
    VecGPU_t vel_gpu(nVec,sh_order);
    { // Copy from host to device
      typedef typename VecGPU_t::device_type device_type;
      trg_gpu.getDevice().Memcpy(trg_gpu.begin(), trg_cpu.begin(),
          trg_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
      src_gpu.getDevice().Memcpy(src_gpu.begin(), src_cpu.begin(),
          src_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
      den_gpu.getDevice().Memcpy(den_gpu.begin(), den_cpu.begin(),
          den_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
      nor_gpu.getDevice().Memcpy(nor_gpu.begin(), nor_cpu.begin(),
          nor_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
      qw_gpu .getDevice().Memcpy( qw_gpu.begin(),  qw_cpu.begin(),
           qw_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
      //vel_gpu.getDevice().Memcpy(vel_gpu.begin(), vel_cpu.begin(),
      //    vel_cpu.size() * sizeof(real), device_type::MemcpyHostToDevice);
    }

    //vel_cpu.getDevice().DirectStokes(src_cpu.begin(), den_cpu.begin(),
    //    qw_cpu.begin(), src_cpu.getStride(), nVec, trg_cpu.begin(), 0, trg_cpu.getStride(), vel_cpu.begin());

    //vel_gpu.getDevice().DirectStokes(src_gpu.begin(), den_gpu.begin(),
    //    qw_gpu.begin(), src_gpu.getStride(), nVec, trg_gpu.begin(), 0, trg_gpu.getStride(), vel_gpu.begin());

    vel_cpu.getDevice().DirectStokesDoubleLayer(src_cpu.begin(), nor_cpu.begin(), den_cpu.begin(),
        qw_cpu.begin(), src_cpu.getStride(), nVec, trg_cpu.begin(), 0, trg_cpu.getStride(), vel_cpu.begin());

    vel_gpu.getDevice().DirectStokesDoubleLayer(src_gpu.begin(), nor_gpu.begin(), den_gpu.begin(),
        qw_gpu.begin(), src_gpu.getStride(), nVec, trg_gpu.begin(), 0, trg_gpu.getStride(), vel_gpu.begin());

    VecCPU_t vel_gpu_copy(nVec,sh_order);
    { // Copy from device to cpu
      typedef typename VecGPU_t::device_type device_type;
      vel_gpu.getDevice().Memcpy(vel_gpu_copy.begin(), vel_gpu.begin(),
          vel_gpu.size() * sizeof(real), device_type::MemcpyDeviceToHost);
    }

    { // Compute error.
      real tol=1e-12;
      VecCPU_t err_vec(nVec,sh_order);
      axpy((real)(-1.0),vel_gpu_copy,vel_cpu,err_vec);
      real err=MaxAbs(err_vec)/MaxAbs(vel_cpu);
      COUT("Error = "<<err<<"\n");
      ASSERT(err<tol,"Velocity error="<<err);
    }
  }

  COUT(emph<<" *** StokesDoubleLayer with GPU device passed ***"<<emph);
#endif //GPU_ACTIVE

  return 0;
}
