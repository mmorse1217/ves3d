#include "PVFMMInterface.h"
#include "Enums.h"

typedef double Real_t;

void nbody(std::vector<Real_t>&  src_coord, std::vector<Real_t>&  src_value, std::vector<Real_t>&  surf_value,
           std::vector<Real_t>&  trg_coord, std::vector<Real_t>&  trg_value,
           const pvfmm::Kernel<double>& kernel_fn, MPI_Comm& comm){
  int np, rank;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &rank);

  long long  n_src =  src_coord.size()/COORD_DIM;
  long long n_trg_glb=0, n_trg = trg_coord.size()/COORD_DIM;
  MPI_Allreduce(&n_trg , & n_trg_glb, 1, MPI_LONG_LONG, MPI_SUM, comm);

  std::vector<Real_t> glb_trg_coord(n_trg_glb*COORD_DIM);
  std::vector<Real_t> glb_trg_value(n_trg_glb*kernel_fn.ker_dim[1],0);
  std::vector<int> recv_disp(np,0);
  { // Gather all target coordinates.
    int send_cnt=n_trg*COORD_DIM;
    std::vector<int> recv_cnts(np);
    MPI_Allgather(&send_cnt    , 1, MPI_INT,
                  &recv_cnts[0], 1, MPI_INT, comm);
    pvfmm::omp_par::scan(&recv_cnts[0], &recv_disp[0], np);
    MPI_Allgatherv(&trg_coord[0]    , send_cnt                    , MPI_DOUBLE,
                   &glb_trg_coord[0], &recv_cnts[0], &recv_disp[0], MPI_DOUBLE, comm);
  }

  { // Evaluate target potential.
    std::vector<Real_t> glb_trg_value_(n_trg_glb*kernel_fn.ker_dim[1],0);
    int omp_p=omp_get_max_threads();
    #pragma omp parallel for
    for(int i=0;i<omp_p;i++){
      size_t a=( i   *n_trg_glb)/omp_p;
      size_t b=((i+1)*n_trg_glb)/omp_p;

      if(kernel_fn.ker_poten!=NULL && src_value.size()>0)
      kernel_fn.ker_poten(&    src_coord[0]            , n_src, &    src_value [0], 1,
                          &glb_trg_coord[0]+a*COORD_DIM,   b-a, &glb_trg_value_[0]+a*kernel_fn.ker_dim[1],NULL);

      if(kernel_fn.dbl_layer_poten!=NULL && surf_value.size()>0)
      kernel_fn.dbl_layer_poten(&    src_coord[0]            , n_src, &   surf_value [0], 1,
                                &glb_trg_coord[0]+a*COORD_DIM,   b-a, &glb_trg_value_[0]+a*kernel_fn.ker_dim[1],NULL);
    }
    MPI_Allreduce(&glb_trg_value_[0], &glb_trg_value[0], glb_trg_value.size(), MPI_DOUBLE, MPI_SUM, comm);
  }

  // Get local target values.
  trg_value.assign(&glb_trg_value[0]+recv_disp[rank]/COORD_DIM*kernel_fn.ker_dim[1], &glb_trg_value[0]+(recv_disp[rank]/COORD_DIM+n_trg)*kernel_fn.ker_dim[1]);
}

int main(int argc, char** argv){
  MPI_Init(&argc,&argv);
  pvfmm::Profile::Enable(true);

  int np, rank;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_size(comm, &np);
  MPI_Comm_rank(comm, &rank);
  srand48(rank);

  size_t N=54076;
  int max_pts=600;
  int mult_order=10;
  int max_depth=20;
  const pvfmm::Kernel<Real_t>* ker=&StokesKernel<Real_t>::Kernel();

  std::vector<Real_t> coord   (N*(COORD_DIM                ));
  std::vector<Real_t> force_sl(N*(    ker->ker_dim[0]));
  std::vector<Real_t> force_dl(N*(COORD_DIM+ker->ker_dim[0]));
  std::vector<Real_t> fmm_v   (N*(    ker->ker_dim[1]));

  for(size_t i=0;i<coord   .size();i++) coord   [i]=drand48();
  for(size_t i=0;i<force_sl.size();i++) force_sl[i]=drand48();
  for(size_t i=0;i<force_dl.size();i++) force_dl[i]=drand48();

  void* ctx=PVFMMCreateContext<Real_t>(-1,max_pts, mult_order, max_depth, ker, comm);

  PVFMMEval(&coord[0], &force_sl[0], &force_dl[0], N, &fmm_v[0], &ctx);
  PVFMMEval(&coord[0], &force_sl[0], &force_dl[0], N, &fmm_v[0], &ctx);

  {// Check error
    size_t n_trg=coord.size()/COORD_DIM;
    size_t N=coord.size()/COORD_DIM*np;

    std::vector<Real_t> trg_sample_coord;
    std::vector<Real_t> trg_sample_value;
    size_t n_trg_sample=0;
    { // Sample target points for verifications.
      size_t n_skip=N*n_trg/1e9;
      if(!n_skip) n_skip=1;
      for(size_t i=0;i<n_trg;i=i+n_skip){
        for(size_t j=0;j<COORD_DIM;j++)
          trg_sample_coord.push_back(coord[i*COORD_DIM+j]);
        for(size_t j=0;j<ker->ker_dim[1];j++)
          trg_sample_value.push_back(fmm_v[i*ker->ker_dim[1]+j]);
        n_trg_sample++;
      }
    }

    // Direct n-body
    std::vector<Real_t> trg_sample_value_(n_trg_sample*ker->ker_dim[1]);
    nbody(           coord, force_sl, force_dl,
          trg_sample_coord, trg_sample_value_, *ker, comm);

    // Compute error
    double max_err=0, max_val=0;
    double max_err_glb=0, max_val_glb=0;
    for(size_t i=0;i<n_trg_sample*ker->ker_dim[1];i++){
      if(fabs(trg_sample_value_[i]-trg_sample_value[i])>max_err)
        max_err=fabs(trg_sample_value_[i]-trg_sample_value[i]);
      if(fabs(trg_sample_value_[i])>max_val)
        max_val=fabs(trg_sample_value_[i]);
    }
    MPI_Reduce(&max_err, &max_err_glb, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&max_val, &max_val_glb, 1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if(!rank) std::cout<<"Maximum Absolute Error:"<<max_err_glb<<'\n';
    if(!rank) std::cout<<"Maximum Relative Error:"<<max_err_glb/max_val_glb<<'\n';
  }

  pvfmm::Profile::print(&comm);
  MPI_Finalize();
  return 0;
}
