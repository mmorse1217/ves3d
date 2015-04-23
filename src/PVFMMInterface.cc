#ifdef HAVE_PVFMM
#include "Enums.h"
#include "Logger.h"

#include <cstring>
#include <parUtils.h>
#include <vector.hpp>
#include <mortonid.hpp>

#include <fmm_pts.hpp>
#include <fmm_node.hpp>
#include <fmm_tree.hpp>

///////////////////////////////////////////////////////////////////////////////
///////////////////////// Kernel Function Declarations ////////////////////////

////////// Stokes Kernel //////////

template <class T>
void stokes_sl_m2l(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(28*dof));
#endif
  const T mu=1.0;
  const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*mu);
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if (R!=0){
        T invR2=1.0/R;
        T invR=sqrt(invR2);
        T invR3=invR2*invR;
        T* f=&v_src[s*4];

        T inner_prod=(f[0]*dR[0] +
                      f[1]*dR[1] +
                      f[2]*dR[2])* invR3;

        T inner_prod_plus_f3_invR3=inner_prod+f[3]*invR3;

        p[0] += f[0]*invR + dR[0]*inner_prod_plus_f3_invR3;
        p[1] += f[1]*invR + dR[1]*inner_prod_plus_f3_invR3;
        p[2] += f[2]*invR + dR[2]*inner_prod_plus_f3_invR3;
      }
    }
    k_out[t*3+0] += p[0]*SCAL_CONST;
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
  }
}

template <class T>
void stokes_m2l_vol_poten(const T* coord, int n, T* out){
  for(int i=0;i<n;i++){
    const T* c=&coord[i*COORD_DIM];
    T rx_2=c[1]*c[1]+c[2]*c[2];
    T ry_2=c[0]*c[0]+c[2]*c[2];
    T rz_2=c[0]*c[0]+c[1]*c[1];
    out[n*3*0+i*3+0]=-rx_2/6; out[n*3*0+i*3+1]=      0; out[n*3*0+i*3+2]=      0;
    out[n*3*1+i*3+0]=      0; out[n*3*1+i*3+1]=-ry_2/6; out[n*3*1+i*3+2]=      0;
    out[n*3*2+i*3+0]=      0; out[n*3*2+i*3+1]=      0; out[n*3*2+i*3+2]=-rz_2/6;
    out[n*3*3+i*3+0]= c[0]/6; out[n*3*3+i*3+1]= c[1]/6; out[n*3*3+i*3+2]= c[2]/6;
  }
}

template <class T>
void stokes_sl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(26*dof));
#endif
  const T mu=1.0;
  const T SCAL_CONST = 1.0/(8.0*pvfmm::const_pi<T>()*mu);
  for(int t=0;t<trg_cnt;t++){
    T p[3]={0,0,0};
    for(int s=0;s<src_cnt;s++){
      T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
               r_trg[3*t+1]-r_src[3*s+1],
               r_trg[3*t+2]-r_src[3*s+2]};
      T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

      T invR2=(R!=0?1.0/R:0.0);
      T invR=sqrt(invR2);
      T invR3=invR2*invR;
      T* f=&v_src[s*3];

      T inner_prod=(f[0]*dR[0] +
                    f[1]*dR[1] +
                    f[2]*dR[2])* invR3;

      p[0] += f[0]*invR + dR[0]*inner_prod;
      p[1] += f[1]*invR + dR[1]*inner_prod;
      p[2] += f[2]*invR + dR[2]*inner_prod;
    }
    k_out[t*3+0] += p[0]*SCAL_CONST;
    k_out[t*3+1] += p[1]*SCAL_CONST;
    k_out[t*3+2] += p[2]*SCAL_CONST;
  }
}

template <class T>
void stokes_dl(T* r_src, int src_cnt, T* v_src, int dof, T* r_trg, int trg_cnt, T* k_out, pvfmm::mem::MemoryManager* mem_mgr){
#ifndef __MIC__
  pvfmm::Profile::Add_FLOP((long long)trg_cnt*(long long)src_cnt*(27*dof));
#endif
  const T SCAL_CONST = -3.0/(4.0*pvfmm::const_pi<T>());
  for(int t=0;t<trg_cnt;t++){
    for(int i=0;i<dof;i++){
      T p[3]={0,0,0};
      for(int s=0;s<src_cnt;s++){
        T dR[3]={r_trg[3*t  ]-r_src[3*s  ],
                 r_trg[3*t+1]-r_src[3*s+1],
                 r_trg[3*t+2]-r_src[3*s+2]};
        T R = (dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);

        if (R!=0){
          T invR2=1.0/R;
          T invR=sqrt(invR2);
          T invR3=invR2*invR;
          T invR5=invR2*invR3;

          T* f=&v_src[(s*dof+i)*6+0];
          T* n=&v_src[(s*dof+i)*6+3];

          T r_dot_n=(n[0]*dR[0]+n[1]*dR[1]+n[2]*dR[2]);
          T r_dot_f=(f[0]*dR[0]+f[1]*dR[1]+f[2]*dR[2]);
          T p_=r_dot_n*r_dot_f*invR5;

          p[0] += dR[0]*p_;
          p[1] += dR[1]*p_;
          p[2] += dR[2]*p_;
        }
      }
      k_out[(t*dof+i)*3+0] += p[0]*SCAL_CONST;
      k_out[(t*dof+i)*3+1] += p[1]*SCAL_CONST;
      k_out[(t*dof+i)*3+2] += p[2]*SCAL_CONST;
    }
  }
}

template <class T>
void stokes_vol_poten(const T* coord, int n, T* out){
  for(int i=0;i<n;i++){
    const T* c=&coord[i*COORD_DIM];
    T rx_2=c[1]*c[1]+c[2]*c[2];
    T ry_2=c[0]*c[0]+c[2]*c[2];
    T rz_2=c[0]*c[0]+c[1]*c[1];
    out[n*3*0+i*3+0]=-rx_2/6; out[n*3*0+i*3+1]=      0; out[n*3*0+i*3+2]=      0;
    out[n*3*1+i*3+0]=      0; out[n*3*1+i*3+1]=-ry_2/6; out[n*3*1+i*3+2]=      0;
    out[n*3*2+i*3+0]=      0; out[n*3*2+i*3+1]=      0; out[n*3*2+i*3+2]=-rz_2/6;
  }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


template<typename T>
void PVFMMBoundingBox(size_t n_src, const T* x, T* scale_xr, T* shift_xr, MPI_Comm comm){
  T& scale_x=*scale_xr;
  T* shift_x= shift_xr;

  if(n_src>0){ // Compute bounding box
    double loc_min_x[DIM];
    double loc_max_x[DIM];
    assert(n_src>0);
    for(size_t k=0;k<DIM;k++){
      loc_min_x[k]=loc_max_x[k]=x[k];
    }

    for(size_t i=0;i<n_src;i++){
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

    T eps=1; // Points should be well within the box.
    while(eps*(T)0.5+(T)1.0>1.0) eps*=0.5; eps*=64;
    scale_x=1.0/(max_x[0]-min_x[0]+2*eps);
    for(size_t k=0;k<DIM;k++){
      scale_x=std::min(scale_x,(T)(1.0/(max_x[k]-min_x[k]+2*eps)));
    }
    if(scale_x*0.0!=0.0) scale_x=1.0; // fix for scal_x=inf
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

  T box_size;
  int max_pts;
  int mult_order;
  int max_depth;
  pvfmm::BoundaryType bndry;
  const pvfmm::Kernel<T>* ker;
  MPI_Comm comm;

  typename Node_t::NodeData tree_data;
  Tree_t* tree;
  Mat_t* mat;
};

template<typename T>
void* PVFMMCreateContext(T box_size, int n, int m, int max_d,
    const pvfmm::Kernel<T>* ker,
    MPI_Comm comm){
  pvfmm::Profile::Tic("FMMContext",&comm,true);
  bool prof_state=pvfmm::Profile::Enable(false);

  // Create new context.
  PVFMMContext<T>* ctx=new PVFMMContext<T>;

  // Set member variables.
  ctx->box_size=box_size;
  ctx->max_pts=n;
  ctx->mult_order=m;
  ctx->max_depth=max_d;
  ctx->bndry=(box_size<=0?pvfmm::FreeSpace:pvfmm::Periodic);
  ctx->ker=ker;
  ctx->comm=comm;

  // Initialize FMM matrices.
  ctx->mat=new typename PVFMMContext<T>::Mat_t();
  ctx->mat->Initialize(ctx->mult_order, ctx->comm, ctx->ker);

  // Set tree_data
  ctx->tree_data.dim=DIM;
  ctx->tree_data.max_depth=ctx->max_depth;
  ctx->tree_data.max_pts=ctx->max_pts;
  { // ctx->tree_data.pt_coord=... //Set points for initial tree.
    int np, myrank;
    MPI_Comm_size(ctx->comm, &np);
    MPI_Comm_rank(ctx->comm, &myrank);

    std::vector<T> coord;
    size_t NN=ceil(pow((T)np*ctx->max_pts,1.0/3.0));
    size_t N_total=NN*NN*NN;
    size_t start= myrank   *N_total/np;
    size_t end  =(myrank+1)*N_total/np;
    for(size_t i=start;i<end;i++){
      coord.push_back(((T)((i/  1    )%NN)+0.5)/NN);
      coord.push_back(((T)((i/ NN    )%NN)+0.5)/NN);
      coord.push_back(((T)((i/(NN*NN))%NN)+0.5)/NN);
    }
    ctx->tree_data.pt_coord=coord;
  }

  // Construct tree.
  bool adap=false; // no data to do adaptive.
  ctx->tree=new typename PVFMMContext<T>::Tree_t(comm);
  ctx->tree->Initialize(&ctx->tree_data);
  ctx->tree->InitFMM_Tree(adap,ctx->bndry);

  pvfmm::Profile::Enable(prof_state);
  pvfmm::Profile::Toc();
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
void PVFMMEval(const T* src_pos, const T* sl_den, size_t n_src, T* trg_vel, void** ctx_){
  PVFMMEval<T>(src_pos, sl_den, NULL, n_src, trg_vel, ctx_);
}

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, T* trg_vel, void** ctx_){
  PVFMMEval(src_pos, sl_den, dl_den, n_src, src_pos, trg_vel, n_src, ctx_);
}

template<typename T>
void PVFMMEval(const T* src_pos, const T* sl_den, const T* dl_den, size_t n_src, const T* trg_pos, T* trg_vel, size_t n_trg, void** ctx_){
  PROFILESTART();
  long long prof_FLOPS=pvfmm::Profile::Add_FLOP(0);
  size_t omp_p=omp_get_max_threads();

  typedef pvfmm::FMM_Node<pvfmm::MPI_Node<T> > Node_t;
  typedef pvfmm::FMM_Pts<Node_t> Mat_t;
  typedef pvfmm::FMM_Tree<Mat_t> Tree_t;

  if(!ctx_[0]) ctx_[0]=PVFMMCreateContext<T>();
  PVFMMContext<T>* ctx=(PVFMMContext<T>*)ctx_[0];
  const int* ker_dim=ctx->ker->ker_dim;

  pvfmm::Profile::Tic("FMM",&ctx->comm);
  T scale_x, shift_x[DIM];
  if(ctx->box_size<=0){ // determine bounding box
    T s0, x0[COORD_DIM];
    T s1, x1[COORD_DIM];
    PVFMMBoundingBox(n_src, src_pos, &s0, x0, ctx->comm);
    PVFMMBoundingBox(n_trg, trg_pos, &s1, x1, ctx->comm);

    T c0[COORD_DIM]={(0.5-x0[0])/s0, (0.5-x0[1])/s0, (0.5-x0[2])/s0};
    T c1[COORD_DIM]={(0.5-x1[0])/s1, (0.5-x1[1])/s1, (0.5-x1[2])/s1};

    scale_x=0;
    scale_x=std::max(scale_x, fabs(c0[0]-c1[0]));
    scale_x=std::max(scale_x, fabs(c0[1]-c1[1]));
    scale_x=std::max(scale_x, fabs(c0[2]-c1[2]));
    scale_x=1.0/(scale_x+1/s0+1/s1);

    shift_x[0]=0.5-(c0[0]+c1[0])*scale_x/2.0;
    shift_x[1]=0.5-(c0[1]+c1[1])*scale_x/2.0;
    shift_x[2]=0.5-(c0[2]+c1[2])*scale_x/2.0;
  }else{
    scale_x=1.0/ctx->box_size;
    shift_x[0]=0;
    shift_x[1]=0;
    shift_x[2]=0;
  }

  pvfmm::Vector<T>  src_scal;
  pvfmm::Vector<T>  trg_scal;
  pvfmm::Vector<T> surf_scal;
  { // Set src_scal, trg_scal
    pvfmm::Vector<T>& src_scal_exp=ctx->ker->src_scal;
    pvfmm::Vector<T>& trg_scal_exp=ctx->ker->trg_scal;
    src_scal .ReInit(ctx->ker->src_scal.Dim());
    trg_scal .ReInit(ctx->ker->trg_scal.Dim());
    surf_scal.ReInit(DIM      +src_scal.Dim());
    for(size_t i=0;i<src_scal.Dim();i++){
      src_scal [i]=pow(scale_x, src_scal_exp[i]);
      surf_scal[i]=scale_x*src_scal[i];
    }
    for(size_t i=0;i<trg_scal.Dim();i++){
      trg_scal[i]=pow(scale_x, trg_scal_exp[i]);
    }
    for(size_t i=src_scal.Dim();i<surf_scal.Dim();i++){
      surf_scal[i]=1;
    }
  }

  pvfmm::Vector<size_t> scatter_index;
  { // Set tree_data
    pvfmm::Vector<T>&  trg_coord=ctx->tree_data. trg_coord;
    pvfmm::Vector<T>&  src_coord=ctx->tree_data. src_coord;
    pvfmm::Vector<T>&  src_value=ctx->tree_data. src_value;
    pvfmm::Vector<T>& surf_value=ctx->tree_data.surf_value;
    pvfmm::Vector<pvfmm::MortonId> pt_mid;

    std::vector<Node_t*> nodes;
    { // Get list of leaf nodes.
      std::vector<Node_t*>& all_nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<all_nodes.size();i++){
        if(all_nodes[i]->IsLeaf() && !all_nodes[i]->IsGhost()){
          nodes.push_back(all_nodes[i]);
        }
      }
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

    { // Set src tree_data
      { // Scatter src data
        // Compute MortonId and copy coordinates and values.
        src_coord .ReInit(       n_src            *DIM   );
        src_value .ReInit(sl_den?n_src*(ker_dim[0]    ):0);
        surf_value.ReInit(dl_den?n_src*(ker_dim[0]+DIM):0);
        pt_mid    .ReInit(n_src);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*n_src)/omp_p;
          size_t b=((tid+1)*n_src)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t j=0;j<DIM;j++){
              src_coord[i*DIM+j]=src_pos[i*DIM+j]*scale_x+shift_x[j];
              while(src_coord[i*DIM+j]< 0.0) src_coord[i*DIM+j]+=1.0;
              while(src_coord[i*DIM+j]>=1.0) src_coord[i*DIM+j]-=1.0;
            }
            pt_mid[i]=pvfmm::MortonId(&src_coord[i*DIM]);
          }
          if(src_value.Dim()) for(size_t i=a;i<b;i++){
            for(size_t j=0;j<ker_dim[0];j++){
              src_value[i*ker_dim[0]+j]=sl_den[i*ker_dim[0]+j]*src_scal[j];
            }
          }
          if(surf_value.Dim()) for(size_t i=a;i<b;i++){
            for(size_t j=0;j<ker_dim[0]+DIM;j++){
              surf_value[i*(ker_dim[0]+DIM)+j]=dl_den[i*(ker_dim[0]+DIM)+j]*surf_scal[j];
            }
          }
        }

        // Scatter src coordinates and values.
        pvfmm::par::SortScatterIndex( pt_mid  , scatter_index, ctx->comm, &min_mid);
        pvfmm::par::ScatterForward  ( pt_mid  , scatter_index, ctx->comm);
        pvfmm::par::ScatterForward  (src_coord, scatter_index, ctx->comm);
        if( src_value.Dim()) pvfmm::par::ScatterForward( src_value, scatter_index, ctx->comm);
        if(surf_value.Dim()) pvfmm::par::ScatterForward(surf_value, scatter_index, ctx->comm);
      }
      { // Set src tree_data
        std::vector<size_t> part_indx(nodes.size()+1);
        part_indx[nodes.size()]=pt_mid.Dim();
        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          part_indx[j]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), nodes[j]->GetMortonId())-&pt_mid[0];
        }

        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          size_t n_pts=part_indx[j+1]-part_indx[j];
          if(src_value.Dim()){
            nodes[j]-> src_coord.ReInit(n_pts*(           DIM),& src_coord[0]+part_indx[j]*(           DIM),false);
            nodes[j]-> src_value.ReInit(n_pts*(ker_dim[0]    ),& src_value[0]+part_indx[j]*(ker_dim[0]    ),false);
          }else{
            nodes[j]-> src_coord.ReInit(0,NULL,false);
            nodes[j]-> src_value.ReInit(0,NULL,false);
          }
          if(surf_value.Dim()){
            nodes[j]->surf_coord.ReInit(n_pts*(           DIM),& src_coord[0]+part_indx[j]*(           DIM),false);
            nodes[j]->surf_value.ReInit(n_pts*(ker_dim[0]+DIM),&surf_value[0]+part_indx[j]*(ker_dim[0]+DIM),false);
          }else{
            nodes[j]->surf_coord.ReInit(0,NULL,false);
            nodes[j]->surf_value.ReInit(0,NULL,false);
          }
        }
      }
    }
    { // Set trg tree_data
      if(trg_pos==src_pos && n_src==n_trg){ // Scatter trg data
        trg_coord.ReInit(src_coord.Dim(),&src_coord[0],false);
      }else{
        // Compute MortonId and copy coordinates.
        trg_coord.Resize(n_trg*DIM);
        pt_mid    .ReInit(n_trg);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*n_trg)/omp_p;
          size_t b=((tid+1)*n_trg)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t j=0;j<DIM;j++){
              trg_coord[i*DIM+j]=trg_pos[i*DIM+j]*scale_x+shift_x[j];
              while(trg_coord[i*DIM+j]< 0.0) trg_coord[i*DIM+j]+=1.0;
              while(trg_coord[i*DIM+j]>=1.0) trg_coord[i*DIM+j]-=1.0;
            }
            pt_mid[i]=pvfmm::MortonId(&trg_coord[i*DIM]);
          }
        }

        // Scatter trg coordinates.
        pvfmm::par::SortScatterIndex( pt_mid  , scatter_index, ctx->comm, &min_mid);
        pvfmm::par::ScatterForward  ( pt_mid  , scatter_index, ctx->comm);
        pvfmm::par::ScatterForward  (trg_coord, scatter_index, ctx->comm);
      }
      { // Set trg tree_data
        std::vector<size_t> part_indx(nodes.size()+1);
        part_indx[nodes.size()]=pt_mid.Dim();
        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          part_indx[j]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), nodes[j]->GetMortonId())-&pt_mid[0];
        }

        #pragma omp parallel for
        for(size_t j=0;j<nodes.size();j++){
          size_t n_pts=part_indx[j+1]-part_indx[j];
          {
            nodes[j]-> trg_coord.ReInit(n_pts*(           DIM),& trg_coord[0]+part_indx[j]*(           DIM),false);
          }
        }
      }
    }
  }

  if(1){ // Optional stuff (redistribute, adaptive refine ...)
    { //Output max tree depth.
      int np, myrank;
      MPI_Comm_size(ctx->comm, &np);
      MPI_Comm_rank(ctx->comm, &myrank);

      long nleaf=0, maxdepth=0;
      std::vector<size_t> all_nodes(MAX_DEPTH+1,0);
      std::vector<size_t> leaf_nodes(MAX_DEPTH+1,0);
      std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<nodes.size();i++){
        Node_t* n=nodes[i];
        if(!n->IsGhost()) all_nodes[n->Depth()]++;
        if(!n->IsGhost() && n->IsLeaf()){
          leaf_nodes[n->Depth()]++;
          if(maxdepth<n->Depth()) maxdepth=n->Depth();
          nleaf++;
        }
      }

      std::stringstream os1,os2;
      os1<<"All Nodes";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=all_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os1<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os1.str());

      os2<<"Leaf Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=leaf_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os2<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os2.str());

      long nleaf_glb=0, maxdepth_glb=0;
      { // MPI_Reduce
        MPI_Allreduce(&nleaf, &nleaf_glb, 1, MPI_INT, MPI_SUM, ctx->comm);
        MPI_Allreduce(&maxdepth, &maxdepth_glb, 1, MPI_INT, MPI_MAX, ctx->comm);
      }
      if(!myrank) COUTDEBUG("Number of Leaf Nodes: "<<nleaf_glb);
      if(!myrank) COUTDEBUG("Tree Depth: "<<maxdepth_glb);
    }
    bool adap=true;
    ctx->tree->InitFMM_Tree(adap,ctx->bndry);
    { //Output max tree depth.
      int np, myrank;
      MPI_Comm_size(ctx->comm, &np);
      MPI_Comm_rank(ctx->comm, &myrank);

      long nleaf=0, maxdepth=0;
      std::vector<size_t> all_nodes(MAX_DEPTH+1,0);
      std::vector<size_t> leaf_nodes(MAX_DEPTH+1,0);
      std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      for(size_t i=0;i<nodes.size();i++){
        Node_t* n=nodes[i];
        if(!n->IsGhost()) all_nodes[n->Depth()]++;
        if(!n->IsGhost() && n->IsLeaf()){
          leaf_nodes[n->Depth()]++;
          if(maxdepth<n->Depth()) maxdepth=n->Depth();
          nleaf++;
        }
      }

      std::stringstream os1,os2;
      os1<<"All  Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=all_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os1<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os1.str());

      os2<<"Leaf Nodes: ";
      for(int i=0;i<MAX_DEPTH;i++){
        int local_size=leaf_nodes[i];
        int global_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_INT, MPI_SUM, ctx->comm);
        os2<<global_size<<' ';
      }
      if(!myrank) COUTDEBUG(os2.str());

      long nleaf_glb=0, maxdepth_glb=0;
      { // MPI_Reduce
        MPI_Allreduce(&nleaf, &nleaf_glb, 1, MPI_INT, MPI_SUM, ctx->comm);
        MPI_Allreduce(&maxdepth, &maxdepth_glb, 1, MPI_INT, MPI_MAX, ctx->comm);
      }
      if(!myrank) COUTDEBUG("Number of Leaf Nodes: "<<nleaf_glb);
      if(!myrank) COUTDEBUG("Tree Depth: "<<maxdepth_glb);
    }
  }

  // Setup tree for FMM.
  ctx->tree->SetupFMM(ctx->mat);
  ctx->tree->RunFMM();

  //Write2File
  //ctx->tree->Write2File("output",0);

  { // Get target potential.
    pvfmm::Vector<T> trg_value;
    { // Get trg data.
      Node_t* n=NULL;
      n=ctx->tree->PreorderFirst();
      while(n!=NULL){
        if(!n->IsGhost() && n->IsLeaf()) break;
        n=ctx->tree->PreorderNxt(n);
      }
      assert(n!=NULL);

      size_t trg_size=0;
      const std::vector<Node_t*>& nodes=ctx->tree->GetNodeList();
      #pragma omp parallel for reduction(+:trg_size)
      for(size_t i=0;i<nodes.size();i++){
        if(nodes[i]->IsLeaf() && !nodes[i]->IsGhost()){
          trg_size+=nodes[i]->trg_value.Dim();
        }
      }
      trg_value.ReInit(trg_size,&n->trg_value[0]);
    }
    pvfmm::par::ScatterReverse  (trg_value, scatter_index, ctx->comm, n_trg);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*n_trg)/omp_p;
      size_t b=((tid+1)*n_trg)/omp_p;
      for(size_t i=a;i<b;i++){
        for(size_t k=0;k<ker_dim[1];k++){
          trg_vel[i*ker_dim[1]+k]=trg_value[i*ker_dim[1]+k]*trg_scal[k];
        }
      }
    }
  }
  pvfmm::Profile::Toc();

  prof_FLOPS=pvfmm::Profile::Add_FLOP(0)-prof_FLOPS;
  //pvfmm::Profile::print(&ctx->comm);
  //PVFMMDestroyContext<T>(ctx_);
  PROFILEEND("",prof_FLOPS);
}

template<typename T>
void PVFMM_GlobalRepart(size_t nv, size_t stride,
    const T* x, const T* tension, size_t* nvr, T** xr,
    T** tensionr, void* user_ptr){
  assert(false); // I think this doesn't work
  exit(0);

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
#endif


