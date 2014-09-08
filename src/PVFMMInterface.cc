#include "Enums.h"

#include <cstring>
#include <parUtils.h>
#include <vector.hpp>
#include <mortonid.hpp>

#include <fmm_pts.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>

template<typename T>
void PVFMMBoundingBox(size_t np, const T* x, T* scale_xr, T* shift_xr, MPI_Comm comm=MPI_COMM_WORLD){
  T& scale_x=*scale_xr;
  T* shift_x= shift_xr;

  if(np>0){ // Compute bounding box
    double loc_min_x[DIM];
    double loc_max_x[DIM];
    assert(np>0);
    for(size_t k=0;k<DIM;k++){
      loc_min_x[k]=loc_max_x[k]=x[k];
    }

    for(size_t i=0;i<np;i++){
      const T* x_=&x[i*DIM];
      for(size_t k=0;k<DIM;k++){
        if(loc_min_x[k]>x_[0]) loc_min_x[k]=x_[0];
        if(loc_max_x[k]<x_[0]) loc_max_x[k]=x_[0];
        ++x_;
      }
    }

    double min_x[DIM];
    double max_x[DIM];
    MPI_Allreduce(loc_min_x, min_x, DIM, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(loc_max_x, max_x, DIM, MPI_DOUBLE, MPI_MAX, comm);

    T eps=0.01; // Points should be well within the box.
    scale_x=1.0/(max_x[0]-min_x[0]);
    for(size_t k=0;k<DIM;k++){
      scale_x=std::min(scale_x,(T)(1.0/(max_x[k]-min_x[k])));
    }
    scale_x*=(1.0-2*eps);
    for(size_t k=0;k<DIM;k++){
      shift_x[k]=-min_x[k]*scale_x+eps;
    }
  }
}

template<typename T>
struct PVFMMContext{
  typedef pvfmm::FMM_Node<pvfmm::MPI_Node<T> > Node_t;
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  int max_pts;
  int mult_order;
  int max_depth;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<T>* ker;
  const pvfmm::Kernel<T>* aux_ker;
  MPI_Comm comm;

  typename Node_t::NodeData tree_data;
  Tree_t* tree;
  Mat_t* mat;
};

template<typename T>
void* PVFMMCreateContext(int n, int m, int max_d,
    pvfmm::BoundaryType bndry,
    const pvfmm::Kernel<T>* ker,
    const pvfmm::Kernel<T>* aux_ker,
    MPI_Comm comm){

  // Create new context.
  PVFMMContext<T>* ctx=new PVFMMContext<T>;

  // Set member variables.
  ctx->max_pts=n;
  ctx->mult_order=m;
  ctx->max_depth=max_d;
  ctx->bndry=bndry;
  ctx->ker=ker;
  ctx->aux_ker=aux_ker;
  ctx->comm=comm;

  // Initialize FMM matrices.
  ctx->mat=new typename PVFMMContext<T>::Mat_t();
  ctx->mat->Initialize(ctx->mult_order, ctx->comm, ctx->ker, ctx->aux_ker);

  // Set tree_data
  ctx->tree_data.dim=DIM;
  ctx->tree_data.max_depth=ctx->max_depth;
  ctx->tree_data.max_pts=ctx->max_pts;
  { // ctx->tree_data.pt_coord=... //Set points for initial tree.
    int np, myrank;
    MPI_Comm_size(ctx->comm, &np);
    MPI_Comm_rank(ctx->comm, &myrank);

    std::vector<double> coord;
    size_t NN=ceil(pow((double)np*ctx->max_pts,1.0/3.0));
    size_t N_total=NN*NN*NN;
    size_t start= myrank   *N_total/np;
    size_t end  =(myrank+1)*N_total/np;
    for(size_t i=start;i<end;i++){
      coord.push_back(((double)((i/  1    )%NN)+0.5)/NN);
      coord.push_back(((double)((i/ NN    )%NN)+0.5)/NN);
      coord.push_back(((double)((i/(NN*NN))%NN)+0.5)/NN);
    }
    ctx->tree_data.pt_coord=coord;
  }

  // Construct tree.
  bool adap=false; // no data to do adaptive.
  ctx->tree=new typename PVFMMContext<T>::Tree_t(comm);
  ctx->tree->Initialize(&ctx->tree_data);
  ctx->tree->InitFMM_Tree(adap,ctx->bndry);

  return ctx;
}

template<typename T>
void PVFMMDestroyContext(void** ctx){
  if(!ctx[0]) return;

  // Delete tree.
  delete ((PVFMMContext<T>*)ctx[0])->tree;

  // Delete matrices.
  delete ((PVFMMContext<T>*)ctx[0])->mat;

  // Delete context.
  delete (PVFMMContext<T>*)ctx[0];
  ctx[0]=NULL;
}

template<typename T>
void PVFMMEval(const T* all_pos, const T* all_den, size_t np, T* all_pot, void** ctx_){
  typedef pvfmm::FMM_Node<pvfmm::MPI_Node<T> > Node_t;
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  if(!ctx_[0]) ctx_[0]=PVFMMCreateContext<T>();
  PVFMMContext<T>* ctx=(PVFMMContext<T>*)ctx_[0];

  T scale_x, shift_x[DIM];
  PVFMMBoundingBox(np, all_pos, &scale_x, shift_x, ctx->comm);

  pvfmm::Vector<size_t> scatter_index;
  pvfmm::Vector<T>& src_coord=ctx->tree_data.src_coord;
  pvfmm::Vector<T>& src_value=ctx->tree_data.src_value;
  { // Set tree_data

    // Compute MortonId and copy coordinates and values.
    src_coord.Resize(np*DIM);
    src_value.Resize(np*DIM);
    pvfmm::Vector<pvfmm::MortonId> src_mid(np);
    #pragma omp parallel for
    for(size_t i=0;i<np;i++){
      for(size_t k=0;k<DIM;k++){
        src_coord[i*DIM+k]=all_pos[i*DIM+k]*scale_x+shift_x[k];
        src_value[i*DIM+k]=all_den[i*DIM+k];
      }
      src_mid[i]=pvfmm::MortonId(&src_coord[i*DIM]);
    }

    pvfmm::MortonId min_mid;
    { // Get first MortonId
      Node_t* n=ctx->tree->PreorderFirst();
      while(n!=NULL){
        if(!n->IsGhost() && n->IsLeaf()) break;
        n=ctx->tree->PreorderNxt(n);
      }
      assert(n!=NULL);
      min_mid=n->GetMortonId();
    }

    // Compute scatter_index.
    pvfmm::par::SortScatterIndex(src_mid  , scatter_index, ctx->comm, &min_mid);

    // Scatter coordinates and values.
    pvfmm::par::ScatterForward  (src_mid  , scatter_index, ctx->comm);
    pvfmm::par::ScatterForward  (src_coord, scatter_index, ctx->comm);
    pvfmm::par::ScatterForward  (src_value, scatter_index, ctx->comm);

    {// Set tree_data
      std::vector<Node_t*> nodes;
      { // Get list of leaf nods.
        std::vector<Node_t*>& all_nodes=ctx->tree->GetNodeList();
        for(size_t i=0;i<all_nodes.size();i++){
          if(all_nodes[i]->IsLeaf() && !all_nodes[i]->IsGhost()){
            nodes.push_back(all_nodes[i]);
          }
        }
      }

      std::vector<size_t> part_indx(nodes.size()+1);
      part_indx[nodes.size()]=src_mid.Dim();
      #pragma omp parallel for
      for(size_t j=0;j<nodes.size();j++){
        part_indx[j]=std::lower_bound(&src_mid[0], &src_mid[0]+src_mid.Dim(), nodes[j]->GetMortonId())-&src_mid[0];
      }

      #pragma omp parallel for
      for(size_t j=0;j<nodes.size();j++){
        size_t n_pts=part_indx[j+1]-part_indx[j];
        nodes[j]->trg_coord.ReInit(n_pts*DIM,&src_coord[0]+part_indx[j]*DIM,false);
        nodes[j]->src_coord.ReInit(n_pts*DIM,&src_coord[0]+part_indx[j]*DIM,false);
        nodes[j]->src_value.ReInit(n_pts*DIM,&src_value[0]+part_indx[j]*DIM,false);
      }
    }
  }

  if(1){ // Optional stuff (redistribute, adaptive refine ...)
    bool adap=true;
    ctx->tree->InitFMM_Tree(adap,ctx->bndry);
  }

  // Setup tree for FMM.
  ctx->tree->SetupFMM(ctx->mat);
  ctx->tree->RunFMM();

  { // Get target potential.
    Node_t* n=NULL;
    { // Get first leaf node.
      n=ctx->tree->PreorderFirst();
      while(n!=NULL){
        if(!n->IsGhost() && n->IsLeaf()) break;
        n=ctx->tree->PreorderNxt(n);
      }
      assert(n!=NULL);
    }
    pvfmm::Vector<T> trg_value(scatter_index.Dim()*DIM,&n->trg_value[0]);

    pvfmm::par::ScatterReverse  (trg_value, scatter_index, ctx->comm, np);
    #pragma omp parallel for
    for(size_t i=0;i<np*DIM;i++){
      all_pot[i]=trg_value[i]*scale_x;
    }
  }
}

template<typename T>
void PVFMM_GlobalRepart(size_t nv, size_t stride,
    const T* x, const T* tension, size_t* nvr, T** xr,
    T** tensionr, void* user_ptr){

  MPI_Comm comm=MPI_COMM_WORLD;

  // Get bounding box.
  T scale_x, shift_x[DIM];
  PVFMMBoundingBox(nv*stride, x, &scale_x, shift_x, comm);

  pvfmm::Vector<pvfmm::MortonId> ves_mid(nv);
  { // Create MortonIds for vesicles.
    scale_x/=stride;
    #pragma omp parallel for
    for(size_t i=0;i<nv;i++){
      T  x_ves[3]={0,0,0};
      const T* x_=&x[i*DIM*stride];
      for(size_t j=0;j<stride;j++){
        for(size_t k=0;k<DIM;k++){
          x_ves[k]+=x_[0];
          ++x_;
        }
      }
      for(size_t k=0;k<DIM;k++){
        x_ves[k]=x_ves[k]*scale_x+shift_x[k];
        assert(x_ves[k]>0.0);
        assert(x_ves[k]<1.0);
      }
      ves_mid[i]=pvfmm::MortonId(x_ves);
    }
  }

  // Determine scatter index vector.
  pvfmm::Vector<size_t> scatter_index;
  pvfmm::par::SortScatterIndex(ves_mid, scatter_index, comm);

  // Allocate memory for output.
  nvr[0]=scatter_index.Dim();
  xr      [0]=new T[nvr[0]*stride*DIM];
  tensionr[0]=new T[nvr[0]*stride    ];

  { // Scatter x
    pvfmm::Vector<T> data(nv*stride*DIM,(T*)x);
    pvfmm::par::ScatterForward(data, scatter_index, comm);

    assert(data.Dim()==nvr[0]*stride*DIM);
    memcpy(xr[0],&data[0],data.Dim()*sizeof(T));
  }

  { // Scatter tension
    pvfmm::Vector<T> data(nv*stride,(T*)tension);
    pvfmm::par::ScatterForward(data, scatter_index, comm);

    assert(data.Dim()==nvr[0]*stride);
    memcpy(tensionr[0],&data[0],data.Dim()*sizeof(T));
  }
}

