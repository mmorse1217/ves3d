#ifdef HAVE_PVFMM

#include <omp.h>
#include <iostream>

#include "DataIO.h"
#include "BgFlow.h"
#include "Array.h"
#include "Vectors.h"
#include "Surface.h"
#include "MovePole.h"
#include "Parameters.h"
#include "OperatorsMats.h"
#include "VesInteraction.h"
#include "PVFMMInterface.h"
#include "InterfacialVelocity.h"

#include <ompUtils.h>
#include <parUtils.h>
#include <mortonid.hpp>
#include <mpi_tree.hpp>
#include <legendre_rule.hpp>
#include <profile.hpp>

template<typename Surf_t>
void NearSingular<Surf_t>::SetSrcCoord(const Surf_t& S_){
  S=&S_;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetTrgCoord(const Surf_t& T_){
  size_t omp_p=omp_get_max_threads();

  const Vec_t& x=T_.getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);

  PVFMMVec_t& pt_coord=T;
  pt_coord.ReInit(N_ves*M_ves*COORD_DIM);

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;

    for(size_t i=a;i<b;i++){
      // read each component of x
      const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
      const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
      const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

      for(size_t j=0;j<M_ves;j++){
        pt_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
        pt_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
        pt_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
      }
    }
  }
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetTrgCoord(const Real_t* trg_coord, size_t N){
  T.ReInit(N*COORD_DIM,trg_coord);
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetSurfaceVel(const Vec_t& S_vel_){
  S_vel=&S_vel_;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetDensity(const Vec_t* force_single_, const Vec_t* force_double_){
  force_single=force_single_;
  force_double=force_double_;
}

template<typename Surf_t>
NearSingular<Surf_t>::Real_t* NearSingular<Surf_t>::EvaluateNearSingular(){
  MPI_Barrier(comm);
  Real_t near=0.5; // TODO: some function of sh_order and accuracy

  int np, rank;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&rank);
  size_t omp_p=omp_get_max_threads();

  struct{
    const Surf_t* S;
    Real_t bbox[4]; // {s,x,y,z} : scale, shift
    Real_t r_near;

    pvfmm::Vector<pvfmm::MortonId> mid; // MortonId of leaf nodes
    pvfmm::Vector<size_t> pt_cnt;       // Point count
    pvfmm::Vector<size_t> pt_dsp;       // Point displ
    pvfmm::MortonId first_mid;          // First non-ghost node

    PVFMMVec_t pt_coord;     // All point coordinates
    pvfmm::Vector<size_t> pt_vesid;     // = pt_id/M_ves
    pvfmm::Vector<size_t> pt_id;        // Scatter id
  } S_let;
  { // Construct S_let
    pvfmm::Profile::Tic("VesLET",&comm,true);
    assert(S); S_let.S=S;

    int np, rank;
    MPI_Comm_size(comm,&np);
    MPI_Comm_rank(comm,&rank);
    size_t omp_p=omp_get_max_threads();

    const Vec_t& x=S->getPosition();
    size_t N_ves = x.getNumSubs(); // Number of vesicles
    size_t M_ves = x.getStride(); // Points per vesicle
    assert(M_ves==x.getGridDim().first*x.getGridDim().second);

    size_t ves_id_offset;
    { // Get ves_id_offset
      long long disp=0;
      long long size=N_ves;
      MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
      ves_id_offset=disp-size;
    }

    size_t tree_depth=0;
    Real_t& r_near=S_let.r_near;
    { // Set tree pt data
      pvfmm::Profile::Tic("PtData",&comm,true);
      Real_t* bbox=S_let.bbox;
      PVFMMVec_t& pt_coord=S_let.pt_coord;
      pvfmm::Vector<size_t>& pt_vesid=S_let.pt_vesid;
      pvfmm::Vector<size_t>& pt_id   =S_let.pt_id   ;

      pt_coord.ReInit(N_ves*M_ves*COORD_DIM);
      pt_vesid.ReInit(N_ves*M_ves);
      pt_id   .ReInit(N_ves*M_ves);

      std::vector<Real_t> r2_ves_(omp_p);
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        size_t a=((tid+0)*N_ves)/omp_p;
        size_t b=((tid+1)*N_ves)/omp_p;

        Real_t r2_ves=0;
        Real_t one_over_M=1.0/M_ves;
        for(size_t i=a;i<b;i++){
          // read each component of x
          const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

          Real_t center_coord[COORD_DIM]={0,0,0};
          for(size_t j=0;j<M_ves;j++){
            center_coord[0]+=xk[j];
            center_coord[1]+=yk[j];
            center_coord[2]+=zk[j];
          }
          center_coord[0]*=one_over_M;
          center_coord[1]*=one_over_M;
          center_coord[2]*=one_over_M;

          for(size_t j=0;j<M_ves;j++){
            pt_coord[(i*M_ves+j)*COORD_DIM+0]=xk[j];
            pt_coord[(i*M_ves+j)*COORD_DIM+1]=yk[j];
            pt_coord[(i*M_ves+j)*COORD_DIM+2]=zk[j];
            pt_vesid[i*M_ves+j]=ves_id_offset+i;

            Real_t dx=(xk[j]-center_coord[0]);
            Real_t dy=(yk[j]-center_coord[1]);
            Real_t dz=(zk[j]-center_coord[2]);
            Real_t r2=dx*dx+dy*dy+dz*dz;
            r2_ves=std::max(r2_ves,r2);
          }
        }
        r2_ves_[tid]=r2_ves;
      }

      Real_t r_ves=0;
      { // Determine r_ves
        double r_ves_loc=0, r_ves_glb=0;
        for(size_t tid=0;tid<omp_p;tid++){
          r_ves_loc=std::max(r2_ves_[tid], r_ves_loc);
        }
        r_ves_loc=sqrt(r_ves_loc);
        MPI_Allreduce(&r_ves_loc, &r_ves_glb, 1, MPI_DOUBLE, MPI_MAX, comm);
        r_ves=r_ves_glb;
      }
      r_near=r_ves*near; // r_near is some function of r_ves.

      { // Determine bbox, tree_depth
        Real_t scale_x, shift_x[COORD_DIM];
        Real_t scale_tmp;
        PVFMMBoundingBox(N_ves*M_ves, &pt_coord[0], &scale_tmp, shift_x, comm);
        { // scale_x, pt_tree_depth
          if(scale_tmp==0){
            scale_tmp=1.0;
            r_near=1.0;
          }
          Real_t domain_length=1.0/scale_tmp+4*r_near;
          Real_t leaf_length=r_near;
          scale_x=1.0/leaf_length;
          while(domain_length*scale_x>1.0){
            scale_x*=0.5;
            tree_depth++;
          }
        }
        for(size_t j=0;j<COORD_DIM;j++){ // Update shift_x
          shift_x[j]=((shift_x[j]/scale_tmp)+2*r_near)*scale_x;
        }
        S_let.bbox[0]=shift_x[0];
        S_let.bbox[1]=shift_x[1];
        S_let.bbox[2]=shift_x[2];
        S_let.bbox[3]=scale_x;
      }
      pvfmm::Profile::Toc();
    }

    pvfmm::Vector<pvfmm::MortonId> mins(np);
    { // Construct local tree
      pvfmm::Profile::Tic("LocalTree",&comm,true);
      pvfmm::Vector<pvfmm::MortonId> pt_mid(N_ves*M_ves);
      PVFMMVec_t& pt_coord=S_let.pt_coord;
      pvfmm::Vector<size_t>& pt_vesid=S_let.pt_vesid;
      pvfmm::Vector<size_t>& pt_id   =S_let.pt_id;

      { // build pt_mid
        Real_t scale_x, shift_x[COORD_DIM];
        { // set scale_x, shift_x
          shift_x[0]=S_let.bbox[0];
          shift_x[1]=S_let.bbox[1];
          shift_x[2]=S_let.bbox[2];
          scale_x=S_let.bbox[3];
        }

        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          Real_t c[3];
          size_t a=((tid+0)*N_ves*M_ves)/omp_p;
          size_t b=((tid+1)*N_ves*M_ves)/omp_p;
          for(size_t i=a;i<b;i++){
            c[0]=pt_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
            c[1]=pt_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
            c[2]=pt_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
            pt_mid[i]=pvfmm::MortonId(c,tree_depth);
          }
        }
      }
      pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm);
      pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
      { // Exchange shared octant with neighbour
        MPI_Allgather(&pt_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                      &  mins[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
        if(rank) assert(mins[rank]!=mins[rank-1]);
        mins[0]=pvfmm::MortonId(0,0,0,tree_depth);
        S_let.first_mid=mins[rank];

        int send_size=0;
        int recv_size=0;
        if(rank<np-1){ // send_size
          send_size=pt_mid.Dim()-(std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), mins[rank+1])-&pt_mid[0]);
        }
        { // recv_size
          MPI_Status status;
          MPI_Sendrecv(&send_size,1,MPI_INT,(rank<np-1?rank+1:0),0,
                       &recv_size,1,MPI_INT,(rank>0?rank-1:np-1),0,comm,&status);
        }

        { // Set new pt_id
          pvfmm::Vector<size_t> pt_id_new(pt_id.Dim()+recv_size-send_size);
          memcpy(&pt_id_new[0]+recv_size, &pt_id[0], (pt_id.Dim()-send_size)*sizeof(size_t));

          MPI_Status status;
          MPI_Sendrecv(&pt_id[0]+pt_id.Dim()-send_size,send_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank<np-1?rank+1:0),0,
                       &pt_id_new[0]                  ,recv_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank>0?rank-1:np-1),0,comm,&status);
          pt_id.Swap(pt_id_new);
        }
        { // Set new pt_mid
          pvfmm::Vector<pvfmm::MortonId> pt_mid_new(pt_mid.Dim()+recv_size-send_size);
          memcpy(&pt_mid_new[0]+recv_size, &pt_mid[0], (pt_mid.Dim()-send_size)*sizeof(pvfmm::MortonId));
          for(size_t i=0;i<recv_size;i++) pt_mid_new[i]=mins[rank];
          pt_mid.Swap(pt_mid_new);
        }
      }
      { // Sort points by pt_id in each octant
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(pt_mid.Dim()*(tid+0))/omp_p;
          size_t b=(pt_mid.Dim()*(tid+1))/omp_p;
          if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
          if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
          for(size_t i=a;i<b;){
            size_t j=i; while(j<b && pt_mid[j]==pt_mid[i]) j++;
            std::sort(&pt_id[0]+i, &pt_id[0]+j);
            i=j;
          }
        }
      }

      pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
      pvfmm::par::ScatterForward(pt_vesid, pt_id, comm);

      pvfmm::Vector<pvfmm::MortonId>& let_mid   =S_let.mid;
      pvfmm::Vector<size_t>&          let_pt_cnt=S_let.pt_cnt;
      pvfmm::Vector<size_t>&          let_pt_dsp=S_let.pt_dsp;
      { // build let_mid, let_pt_cnt, let_pt_dsp
        std::vector<std::vector<pvfmm::MortonId> > mid_(omp_p);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          std::vector<pvfmm::MortonId>& mid=mid_[tid];
          size_t a=(pt_mid.Dim()*(tid+0))/omp_p;
          size_t b=(pt_mid.Dim()*(tid+1))/omp_p;
          if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
          if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
          if(a<b) mid.push_back(pt_mid[a]);
          for(size_t i=a;i<b;i++){
            if(mid.back()!=pt_mid[i]) mid.push_back(pt_mid[i]);
          }
        }
        { // Resize let_mid
          size_t size=0;
          for(size_t tid=0;tid<omp_p;tid++){
            size+=mid_[tid].size();
          }
          let_mid.ReInit(size);
        }
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){ // Set let_mid
          size_t offset=0;
          for(size_t i=0;i<tid;i++){
            offset+=mid_[i].size();
          }

          std::vector<pvfmm::MortonId>& mid=mid_[tid];
          for(size_t i=0;i<mid.size();i++){
            let_mid[i+offset]=mid[i];
          }
        }

        // Resize let_pt_cnt, let_pt_dsp
        let_pt_cnt.ReInit(let_mid.Dim());
        let_pt_dsp.ReInit(let_mid.Dim());

        #pragma omp parallel for
        for(size_t i=0;i<let_mid.Dim();i++){
          let_pt_dsp[i]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mid[i])-&pt_mid[0];
        }

        #pragma omp parallel for
        for(size_t i=1;i<let_mid.Dim();i++){
          let_pt_cnt[i-1]=let_pt_dsp[i]-let_pt_dsp[i-1];
        }
        if(let_mid.Dim())
        let_pt_cnt[let_mid.Dim()-1]=pt_mid.Dim()-let_pt_dsp[let_mid.Dim()-1];
      }
      pvfmm::Profile::Toc();
    }

    { // Construct LET
      pvfmm::Profile::Tic("LETree",&comm,true);
      // TODO: Replace MPI_Alltoallv with Mpi_Alltoallv_sparse

      pvfmm::Vector<size_t> snode_id;
      pvfmm::Vector<int> snode_cnt(np);
      pvfmm::Vector<int> snode_dsp(np);
      pvfmm::Vector<int> rnode_cnt(np);
      pvfmm::Vector<int> rnode_dsp(np);
      { // Compute shared_pair
        pvfmm::Vector<pvfmm::par::SortPair<int,size_t> > shared_pair; // pid, node_id list
        { // Compute shared_pair
          std::vector<std::vector<pvfmm::par::SortPair<int,size_t> > > shared_pair_(omp_p);
          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){
            pvfmm::Vector<pvfmm::MortonId>& let_mid=S_let.mid;

            Real_t coord[3];
            Real_t s=pow(0.5,tree_depth);
            pvfmm::par::SortPair<int, size_t> pair;
            size_t a=(let_mid.Dim()*(tid+0))/omp_p;
            size_t b=(let_mid.Dim()*(tid+1))/omp_p;
            std::set<int> pid_set;
            for(size_t i=a;i<b;i++){
              pid_set.clear();
              let_mid[i].GetCoord(coord);
              for(int j0=-1;j0<=1;j0++)
              for(int j1=-1;j1<=1;j1++)
              for(int j2=-1;j2<=1;j2++)
              if(j0 || j1 || j2){
                Real_t c[3]={coord[0]+j0*s+0.5*s,coord[1]+j1*s+0.5*s,coord[2]+j2*s+0.5*s};
                if(c[0]>0 && c[1]>0 && c[2]>0 && c[0]<1.0 && c[1]<1.0 && c[2]<1.0){
                  pvfmm::MortonId m(c,tree_depth+1);
                  if((rank && m<mins[rank]) || (rank<np-1 && m>=mins[rank+1])){
                    int pid=std::lower_bound(&mins[0], &mins[0]+np, m)-&mins[0]-1;
                    assert(pid!=rank);
                    assert(pid>=0);
                    assert(pid<np);
                    pid_set.insert(pid);
                  }
                }
              }
              for(std::set<int>::iterator it=pid_set.begin(); it!=pid_set.end(); ++it){ // Add shared pair
                pair.data=i;
                pair.key=*it;
                shared_pair_[tid].push_back(pair);
              }
            }
          }
          { // Resize shared_pair
            size_t size=0;
            for(size_t tid=0;tid<omp_p;tid++){
              size+=shared_pair_[tid].size();
            }
            shared_pair.ReInit(size);
          }
          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){
            size_t offset=0;
            for(size_t i=0;i<tid;i++){
              offset+=shared_pair_[i].size();
            }
            for(size_t i=0;i<shared_pair_[tid].size();i++){
              shared_pair[offset+i]=shared_pair_[tid][i];
            }
          }
          pvfmm::omp_par::merge_sort(&shared_pair[0],&shared_pair[0]+shared_pair.Dim());
        }
        { // Set snode_cnt, snode_idx
          snode_cnt.SetZero();
          snode_id.ReInit(shared_pair.Dim());
          for(size_t i=0;i<shared_pair.Dim();i++){
            snode_cnt[shared_pair[i].key]++;
            snode_id[i]=shared_pair[i].data;
          }
        }

        // Get rnode_cnt
        MPI_Alltoall(&snode_cnt[0], 1, pvfmm::par::Mpi_datatype<int>::value(),
                     &rnode_cnt[0], 1, pvfmm::par::Mpi_datatype<int>::value(), comm);

        snode_dsp[0]=0; pvfmm::omp_par::scan(&snode_cnt[0], &snode_dsp[0], snode_cnt.Dim());
        rnode_dsp[0]=0; pvfmm::omp_par::scan(&rnode_cnt[0], &rnode_dsp[0], rnode_cnt.Dim());
      }


      pvfmm::Vector<pvfmm::MortonId> recv_mid;
      pvfmm::Vector<size_t> recv_pt_cnt;
      pvfmm::Vector<size_t> recv_pt_dsp;

      pvfmm::Vector<pvfmm::MortonId> send_mid;
      pvfmm::Vector<size_t> send_pt_cnt;
      pvfmm::Vector<size_t> send_pt_dsp;
      { // Send-recv node data
        size_t send_size=snode_cnt[np-1]+snode_dsp[np-1];
        send_mid.ReInit(send_size);
        send_pt_cnt.ReInit(send_size);
        send_pt_dsp.ReInit(send_size+1);

        size_t recv_size=rnode_cnt[np-1]+rnode_dsp[np-1];
        recv_mid.ReInit(recv_size);
        recv_pt_cnt.ReInit(recv_size);
        recv_pt_dsp.ReInit(recv_size+1);

        for(size_t i=0;i<send_size;i++){ // Set send data
          send_mid   [i]=S_let.mid   [snode_id[i]];
          send_pt_cnt[i]=S_let.pt_cnt[snode_id[i]];
        }
        MPI_Alltoallv(&send_mid[0], &snode_cnt[0], &snode_dsp[0], pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                      &recv_mid[0], &rnode_cnt[0], &rnode_dsp[0], pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
        MPI_Alltoallv(&send_pt_cnt[0], &snode_cnt[0], &snode_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(),
                      &recv_pt_cnt[0], &rnode_cnt[0], &rnode_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(), comm);

        send_pt_dsp[0]=0; pvfmm::omp_par::scan(&send_pt_cnt[0], &send_pt_dsp[0], send_pt_cnt.Dim()+1);
        recv_pt_dsp[0]=0; pvfmm::omp_par::scan(&recv_pt_cnt[0], &recv_pt_dsp[0], recv_pt_cnt.Dim()+1);
      }


      PVFMMVec_t recv_pt_coord;
      pvfmm::Vector<size_t> recv_pt_vesid;
      pvfmm::Vector<size_t> recv_pt_id;

      PVFMMVec_t send_pt_coord;
      pvfmm::Vector<size_t> send_pt_vesid;
      pvfmm::Vector<size_t> send_pt_id;
      { // Send-recv pt data
        size_t send_size_pt=send_pt_dsp[send_mid.Dim()];
        send_pt_coord.ReInit(send_size_pt*COORD_DIM);
        send_pt_vesid.ReInit(send_size_pt);
        send_pt_id   .ReInit(send_size_pt);
        { // Set send data
          #pragma omp parallel for
          for(size_t i=0;i<send_pt_cnt.Dim();i++){
            size_t offset_in=S_let.pt_dsp[snode_id[i]];
            size_t offset_out=send_pt_dsp[i];
            for(size_t j=0;j<send_pt_cnt[i];j++){
              send_pt_coord[(offset_out+j)*COORD_DIM+0]=S_let.pt_coord[(offset_in+j)*COORD_DIM+0];
              send_pt_coord[(offset_out+j)*COORD_DIM+1]=S_let.pt_coord[(offset_in+j)*COORD_DIM+1];
              send_pt_coord[(offset_out+j)*COORD_DIM+2]=S_let.pt_coord[(offset_in+j)*COORD_DIM+2];
              send_pt_vesid[ offset_out+j             ]=S_let.pt_vesid[ offset_in+j             ];
              send_pt_id   [ offset_out+j             ]=S_let.pt_id   [ offset_in+j             ];
            }
          }
        }

        size_t recv_size_pt=recv_pt_dsp[recv_mid.Dim()];
        recv_pt_coord.ReInit(recv_size_pt*COORD_DIM);
        recv_pt_vesid.ReInit(recv_size_pt);
        recv_pt_id   .ReInit(recv_size_pt);

        { // Send-recv data
          pvfmm::Vector<int> send_cnt(np), send_dsp(np);
          pvfmm::Vector<int> recv_cnt(np), recv_dsp(np);
          { // Set send_dsp, send_cnt,  recv_dsp, recv_cnt
            #pragma omp parallel for
            for(size_t i=0;i<np;i++){
              send_dsp[i]=send_pt_dsp[snode_dsp[i]];
              recv_dsp[i]=recv_pt_dsp[rnode_dsp[i]];
            }
            #pragma omp parallel for
            for(size_t i=0;i<np-1;i++){
              send_cnt[i]=send_dsp[i+1]-send_dsp[i];
              recv_cnt[i]=recv_dsp[i+1]-recv_dsp[i];
            }
            send_cnt[np-1]=send_size_pt-send_dsp[np-1];
            recv_cnt[np-1]=recv_size_pt-recv_dsp[np-1];
          }

          MPI_Alltoallv(&send_pt_id   [0], &send_cnt[0], &send_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(),
                        &recv_pt_id   [0], &recv_cnt[0], &recv_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(), comm);
          MPI_Alltoallv(&send_pt_vesid[0], &send_cnt[0], &send_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(),
                        &recv_pt_vesid[0], &recv_cnt[0], &recv_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(), comm);
          for(size_t i=0;i<np;i++){
            send_cnt[i]*=COORD_DIM; send_dsp[i]*=COORD_DIM;
            recv_cnt[i]*=COORD_DIM; recv_dsp[i]*=COORD_DIM;
          }
          MPI_Alltoallv(&send_pt_coord[0], &send_cnt[0], &send_dsp[0], pvfmm::par::Mpi_datatype<Real_t>::value(),
                        &recv_pt_coord[0], &recv_cnt[0], &recv_dsp[0], pvfmm::par::Mpi_datatype<Real_t>::value(), comm);
        }
      }

      { // Add ghost nodes to S_let
        pvfmm::Vector<pvfmm::MortonId> new_mid(S_let.mid .Dim()+recv_mid     .Dim());
        pvfmm::Vector<size_t> new_pt_cnt  (S_let.pt_cnt  .Dim()+recv_pt_cnt  .Dim());
        pvfmm::Vector<size_t> new_pt_dsp  (S_let.pt_dsp  .Dim()+recv_pt_dsp  .Dim());

        PVFMMVec_t new_pt_coord(S_let.pt_coord.Dim()+recv_pt_coord.Dim());
        pvfmm::Vector<size_t> new_pt_vesid(S_let.pt_vesid.Dim()+recv_pt_vesid.Dim());
        pvfmm::Vector<size_t> new_pt_id   (S_let.pt_id   .Dim()+recv_pt_id   .Dim());

        { // Copy mid
          size_t offset=0, size=0;
          size_t recv_split=rnode_dsp[rank];
          size=                recv_split; memcpy(&new_mid[0]+offset, & recv_mid[0]           , size*sizeof(pvfmm::MortonId)); offset+=size;
          size=S_let.mid.Dim()           ; memcpy(&new_mid[0]+offset, &S_let.mid[0]           , size*sizeof(pvfmm::MortonId)); offset+=size;
          size= recv_mid.Dim()-recv_split; memcpy(&new_mid[0]+offset, & recv_mid[0]+recv_split, size*sizeof(pvfmm::MortonId));
        }
        { // Copy pt_cnt
          size_t offset=0, size=0;
          size_t recv_split=rnode_dsp[rank];
          size=                   recv_split; memcpy(&new_pt_cnt[0]+offset, & recv_pt_cnt[0]           , size*sizeof(size_t)); offset+=size;
          size=S_let.pt_cnt.Dim()           ; memcpy(&new_pt_cnt[0]+offset, &S_let.pt_cnt[0]           , size*sizeof(size_t)); offset+=size;
          size= recv_pt_cnt.Dim()-recv_split; memcpy(&new_pt_cnt[0]+offset, & recv_pt_cnt[0]+recv_split, size*sizeof(size_t));
        }
        { // Compute pt_dsp
          new_pt_dsp[0]=0; pvfmm::omp_par::scan(&new_pt_cnt[0], &new_pt_dsp[0], new_pt_cnt.Dim());
        }

        { // Copy pt_coord
          size_t offset=0, size=0;
          size_t recv_split=recv_pt_dsp[rnode_dsp[rank]]*COORD_DIM;
          size=                     recv_split; memcpy(&new_pt_coord[0]+offset, & recv_pt_coord[0]           , size*sizeof(Real_t)); offset+=size;
          size=S_let.pt_coord.Dim()           ; memcpy(&new_pt_coord[0]+offset, &S_let.pt_coord[0]           , size*sizeof(Real_t)); offset+=size;
          size= recv_pt_coord.Dim()-recv_split; memcpy(&new_pt_coord[0]+offset, & recv_pt_coord[0]+recv_split, size*sizeof(Real_t));
        }
        { // Copy pt_vesid
          size_t offset=0, size=0;
          size_t recv_split=recv_pt_dsp[rnode_dsp[rank]];
          size=                     recv_split; memcpy(&new_pt_vesid[0]+offset, & recv_pt_vesid[0]           , size*sizeof(size_t)); offset+=size;
          size=S_let.pt_vesid.Dim()           ; memcpy(&new_pt_vesid[0]+offset, &S_let.pt_vesid[0]           , size*sizeof(size_t)); offset+=size;
          size= recv_pt_vesid.Dim()-recv_split; memcpy(&new_pt_vesid[0]+offset, & recv_pt_vesid[0]+recv_split, size*sizeof(size_t));
        }
        { // Copy pt_id
          size_t offset=0, size=0;
          size_t recv_split=recv_pt_dsp[rnode_dsp[rank]];
          size=                  recv_split; memcpy(&new_pt_id[0]+offset, & recv_pt_id[0]           , size*sizeof(size_t)); offset+=size;
          size=S_let.pt_id.Dim()           ; memcpy(&new_pt_id[0]+offset, &S_let.pt_id[0]           , size*sizeof(size_t)); offset+=size;
          size= recv_pt_id.Dim()-recv_split; memcpy(&new_pt_id[0]+offset, & recv_pt_id[0]+recv_split, size*sizeof(size_t));
        }

        new_mid     .Swap(S_let.mid     );
        new_pt_cnt  .Swap(S_let.pt_cnt  );
        new_pt_dsp  .Swap(S_let.pt_dsp  );

        new_pt_coord.Swap(S_let.pt_coord);
        new_pt_vesid.Swap(S_let.pt_vesid);
        new_pt_id   .Swap(S_let.pt_id   );
      }
      pvfmm::Profile::Toc();
    }

    pvfmm::Profile::Toc();
  }

  PVFMMVec_t near_trg_coord;
  pvfmm::Vector<size_t> near_trg_cnt;
  pvfmm::Vector<size_t> trg_scatter; // Scatter velocity to original location
  pvfmm::Vector<size_t> trg_pt_id;   // target index at original location
  { // find vesicle, near target pairs
    pvfmm::Profile::Tic("SetTrgCoord",&comm,true);
    pvfmm::Vector<size_t> pt_id;
    PVFMMVec_t pt_coord=T;
    size_t trg_id_offset;
    { // Get trg_id_offset
      long long disp=0;
      long long size=pt_coord.Dim()/COORD_DIM;
      MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
      trg_id_offset=disp-size;
    }

    pvfmm::Vector<pvfmm::MortonId> tree_mid;
    pvfmm::Vector<size_t>          tree_pt_cnt;
    pvfmm::Vector<size_t>          tree_pt_dsp;
    { // Construct target point tree
      pvfmm::Profile::Tic("TrgTree",&comm,true);
      pvfmm::Vector<pvfmm::MortonId> pt_mid(pt_coord.Dim()/COORD_DIM);

      { // Set pt_mid
        Real_t scale_x, shift_x[COORD_DIM];
        { // set scale_x, shift_x
          shift_x[0]=S_let.bbox[0];
          shift_x[1]=S_let.bbox[1];
          shift_x[2]=S_let.bbox[2];
          scale_x=S_let.bbox[3];
        }
        assert(S_let.mid.Dim());
        size_t tree_depth=S_let.mid[0].GetDepth();
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          Real_t c[3];
          size_t a=((tid+0)*T.Dim()/COORD_DIM)/omp_p;
          size_t b=((tid+1)*T.Dim()/COORD_DIM)/omp_p;
          for(size_t i=a;i<b;i++){
            c[0]=pt_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
            c[1]=pt_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
            c[2]=pt_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
            assert(c[0]>0.0 && c[0]<1.0);
            assert(c[1]>0.0 && c[1]<1.0);
            assert(c[2]>0.0 && c[2]<1.0);
            pt_mid[i]=pvfmm::MortonId(c,tree_depth);
          }
        }
      }

      { // Sort pt data (pt_mid, pt_coord, pt_id)
        pvfmm::MortonId first_mid=S_let.first_mid;
        pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm, &first_mid);
        pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
        pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
      }

      { // build tree_mid, tree_pt_cnt, tree_pt_dsp
        std::vector<std::vector<pvfmm::MortonId> > mid_(omp_p);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          std::vector<pvfmm::MortonId>& mid=mid_[tid];
          size_t a=(pt_mid.Dim()*(tid+0))/omp_p;
          size_t b=(pt_mid.Dim()*(tid+1))/omp_p;
          if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
          if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
          if(a<b) mid.push_back(pt_mid[a]);
          for(size_t i=a;i<b;i++){
            if(mid.back()!=pt_mid[i]) mid.push_back(pt_mid[i]);
          }
        }
        { // Resize tree_mid
          size_t size=0;
          for(size_t tid=0;tid<omp_p;tid++){
            size+=mid_[tid].size();
          }
          tree_mid.ReInit(size);
        }
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){ // Set tree_mid
          size_t offset=0;
          for(size_t i=0;i<tid;i++){
            offset+=mid_[i].size();
          }

          std::vector<pvfmm::MortonId>& mid=mid_[tid];
          for(size_t i=0;i<mid.size();i++){
            tree_mid[i+offset]=mid[i];
          }
        }

        // Resize tree_pt_cnt, tree_pt_dsp
        tree_pt_cnt.ReInit(tree_mid.Dim());
        tree_pt_dsp.ReInit(tree_mid.Dim());

        #pragma omp parallel for
        for(size_t i=0;i<tree_mid.Dim();i++){
          tree_pt_dsp[i]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), tree_mid[i])-&pt_mid[0];
        }

        #pragma omp parallel for
        for(size_t i=1;i<tree_mid.Dim();i++){
          tree_pt_cnt[i-1]=tree_pt_dsp[i]-tree_pt_dsp[i-1];
        }
        if(tree_mid.Dim())
        tree_pt_cnt[tree_mid.Dim()-1]=pt_mid.Dim()-tree_pt_dsp[tree_mid.Dim()-1];
      }
      pvfmm::Profile::Toc();
    }

    pvfmm::Vector<size_t> near_trg_pt_id;
    pvfmm::Vector<size_t> near_ves_pt_id;
    // near_ves_pt_id: The neares vesicle point to the target point
    // TODO: This is not used right now, and the neares vesicle point is
    // compute again when constructing quadratic surface patch

    { // Find near trg, ves points
      pvfmm::Profile::Tic("NearPair",&comm,true);
      std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_(omp_p); // (loc_trg_id, S_let.pt_id)
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_[tid];
        size_t tree_depth; Real_t r2_near;
        { // Set tree_depth, r_near
          tree_depth=S_let.mid[0].GetDepth();
          r2_near=S_let.r_near;
          r2_near*=r2_near;
        }
        Real_t s=std::pow(0.5,tree_depth);

        size_t FLOP=0;
        size_t a=((tid+0)*tree_mid.Dim())/omp_p;
        size_t b=((tid+1)*tree_mid.Dim())/omp_p;
        for(size_t i=a;i<b;i++){
          size_t tcnt=tree_pt_cnt[i];
          size_t tdsp=tree_pt_dsp[i];
          PVFMMVec_t tcoord;
          { // Set t_coord
            tcoord.ReInit(tcnt*COORD_DIM,&pt_coord[tdsp*COORD_DIM],false);
          }

          PVFMMVec_t scoord[27];
          pvfmm::Vector<size_t> spt_id[27];
          pvfmm::Vector<size_t> svesid[27];
          { // Set scoord, spt_id
            for(size_t j=0;j<27;j++){
              scoord[j].ReInit(0);
              spt_id[j].ReInit(0);
            }
            size_t indx=0;
            Real_t coord[COORD_DIM];
            tree_mid[i].GetCoord(coord);
            for(int j0=-1;j0<=1;j0++)
            for(int j1=-1;j1<=1;j1++)
            for(int j2=-1;j2<=1;j2++){
              Real_t c[3]={coord[0]+j0*s+0.5*s,coord[1]+j1*s+0.5*s,coord[2]+j2*s+0.5*s};
              if(c[0]>0 && c[1]>0 && c[2]>0 && c[0]<1.0 && c[1]<1.0 && c[2]<1.0){
                pvfmm::MortonId m(c,tree_depth);
                pvfmm::Vector<pvfmm::MortonId>& mid=S_let.mid;
                int k=std::lower_bound(&mid[0], &mid[0]+mid.Dim(), m)-&mid[0];
                if(k<mid.Dim() && mid[k]==m){
                  scoord[indx].ReInit(S_let.pt_cnt[k]*COORD_DIM, &S_let.pt_coord[0]+S_let.pt_dsp[k]*COORD_DIM, false);
                  spt_id[indx].ReInit(S_let.pt_cnt[k]          , &S_let.pt_id   [0]+S_let.pt_dsp[k]          , false);
                  svesid[indx].ReInit(S_let.pt_cnt[k]          , &S_let.pt_vesid[0]+S_let.pt_dsp[k]          , false);
                }
              }
              indx++;
            }
          }

          { // Find near pairs
            pvfmm::par::SortPair<size_t, std::pair<size_t, Real_t> > pair_data;
            std::vector<pvfmm::par::SortPair<size_t, std::pair<size_t, Real_t> > > pair_vec;
            for(size_t t=0;t<tcnt;t++){
              pair_vec.clear();
              for(size_t k=0;k<27;k++)
              for(size_t s=0;s<spt_id[k].Dim();s++){
                Real_t dx=scoord[k][s*COORD_DIM+0]-tcoord[t*COORD_DIM+0];
                Real_t dy=scoord[k][s*COORD_DIM+1]-tcoord[t*COORD_DIM+1];
                Real_t dz=scoord[k][s*COORD_DIM+2]-tcoord[t*COORD_DIM+2];
                Real_t r2=dx*dx+dy*dy+dz*dz;
                if(r2<r2_near){
                  size_t vesid=svesid[k][s];
                  size_t pt_id=spt_id[k][s];

                  pair_data.key        =vesid;
                  pair_data.data.first =pt_id;
                  pair_data.data.second=r2;
                  if(pair_vec.size() && pair_vec.back().key==vesid){
                    std::pair<size_t, Real_t>& back_data=pair_vec.back().data;
                    if(r2<back_data.second) back_data=pair_data.data;
                  }else{
                    pair_vec.push_back(pair_data);
                  }
                }
              }
              std::sort(&pair_vec[0],&pair_vec[0]+pair_vec.size());

              Real_t r2min=0;
              for(size_t s=0;s<pair_vec.size();s++){
                std::pair<size_t, size_t> new_pair;
                new_pair.first=tdsp+t;
                new_pair.second=pair_vec[s].data.first;
                if(s && pair_vec[s-1].key==pair_vec[s].key){
                  Real_t r2=pair_vec[s].data.second;
                  if(r2<r2min){
                    r2min=r2;
                    near_pair.back()=new_pair;
                  }
                }else{
                  r2min=pair_vec[s].data.second;
                  near_pair.push_back(new_pair);
                }
              }
            }
            for(size_t k=0;k<27;k++){ // add FLOPS
              FLOP+=spt_id[k].Dim()*tcnt*8;
            }
          }
        }
        pvfmm::Profile::Add_FLOP(FLOP);
      }

      size_t near_size=0;
      pvfmm::Vector<size_t> near_cnt(omp_p);
      pvfmm::Vector<size_t> near_dsp(omp_p); near_dsp[0]=0;
      for(size_t tid=0;tid<omp_p;tid++){
        if(tid)
        near_dsp[tid]=near_pair_[tid-1].size()+near_dsp[tid-1];
        near_cnt[tid]=near_pair_[tid  ].size();
        near_size   +=near_pair_[tid  ].size();
      }

      near_trg_coord.ReInit(near_size*COORD_DIM);
      near_trg_pt_id.ReInit(near_size          );
      near_ves_pt_id.ReInit(near_size          );

      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        size_t dsp=near_dsp[tid];
        size_t cnt=near_cnt[tid];
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_[tid];
        for(size_t i=0;i<cnt;i++){
          size_t loc_trg_id=near_pair[i].first;
          near_trg_coord[(dsp+i)*COORD_DIM+0]=pt_coord[loc_trg_id*COORD_DIM+0];
          near_trg_coord[(dsp+i)*COORD_DIM+1]=pt_coord[loc_trg_id*COORD_DIM+1];
          near_trg_coord[(dsp+i)*COORD_DIM+2]=pt_coord[loc_trg_id*COORD_DIM+2];
          near_trg_pt_id[(dsp+i)            ]=pt_id[loc_trg_id];
          near_ves_pt_id[(dsp+i)            ]=near_pair[i].second;
        }
      }
      pvfmm::Profile::Toc();
    }

    { // Scatter trg points to vesicles
      pvfmm::Profile::Tic("ScatterTrg",&comm,true);
      const Vec_t& x=S_let.S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);

      size_t ves_pt_id_offset;
      { // Get ves_pt_id_offset
        long long disp=0;
        long long size=N_ves;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        ves_pt_id_offset=(disp-size)*M_ves;
      }

      pvfmm::Vector<size_t> scatter_idx;
      pvfmm::par::SortScatterIndex(near_ves_pt_id, scatter_idx, comm,&ves_pt_id_offset);
      pvfmm::par::ScatterForward  (near_ves_pt_id, scatter_idx, comm);
      pvfmm::par::ScatterForward  (near_trg_coord, scatter_idx, comm);
      pvfmm::par::ScatterForward  (near_trg_pt_id, scatter_idx, comm);

      { // Compute near_trg_cnt
        near_trg_cnt.ReInit(N_ves);near_trg_cnt.SetZero();
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(near_trg_pt_id.Dim()*(tid+0))/omp_p;
          size_t b=(near_trg_pt_id.Dim()*(tid+1))/omp_p;
          while(a>0 && a<near_trg_pt_id.Dim() && near_ves_pt_id[a-1]/M_ves==near_ves_pt_id[a]/M_ves) a++;
          while(b>0 && b<near_trg_pt_id.Dim() && near_ves_pt_id[b-1]/M_ves==near_ves_pt_id[b]/M_ves) b++;
          for(size_t i=a;i<b;i++) near_trg_cnt[(near_ves_pt_id[i]-ves_pt_id_offset)/M_ves]++;
        }
      }
      pvfmm::Profile::Toc();
    }

    { // Determine trg_scatter, trg_pt_id
      trg_pt_id.Swap(near_trg_pt_id);
      pvfmm::par::SortScatterIndex(trg_pt_id, trg_scatter, comm, &trg_id_offset);
      pvfmm::par::ScatterForward  (trg_pt_id, trg_scatter, comm);
    }
    pvfmm::Profile::Toc();
  }

  { // Compute near-singular
    pvfmm::Profile::Tic("NearSing",&comm,true);
    pvfmm::Vector<size_t> near_trg_dsp;
    { // Compute near_trg_dsp
      size_t N_ves=near_trg_cnt.Dim();
      near_trg_dsp.ReInit(N_ves);near_trg_dsp[0]=0;
      pvfmm::omp_par::scan(&near_trg_cnt[0], &near_trg_dsp[0], N_ves);
      assert(near_trg_dsp[N_ves-1]+near_trg_cnt[N_ves-1]==near_trg_coord.Dim()/COORD_DIM);
    }
    trg_vel.ReInit(near_trg_coord.Dim()); trg_vel.SetZero();
    StokesNearSingular(S_let.r_near, &near_trg_cnt[0], &near_trg_dsp[0], &near_trg_coord[0], &trg_vel[0]);
    pvfmm::Profile::Toc();
  }

  if(trg_pt_id.Dim()){ // Scatter trg velocity
    pvfmm::Profile::Tic("ScatterTrg",&comm,true);
    size_t trg_id_offset=trg_pt_id[0];
    size_t trg_count=trg_pt_id[trg_pt_id.Dim()-1]-trg_id_offset+1;
    assert(T.Dim()==trg_count*COORD_DIM);

    pvfmm::par::ScatterForward(trg_vel, trg_scatter, comm);
    assert(trg_vel.Dim()==trg_pt_id.Dim()*COORD_DIM);

    PVFMMVec_t trg_vel_final(trg_count*COORD_DIM);
    trg_vel_final.SetZero();

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*trg_pt_id.Dim())/omp_p;
      size_t b=((tid+1)*trg_pt_id.Dim())/omp_p;
      while(a>0 && a<trg_pt_id.Dim() && trg_pt_id[a-1]==trg_pt_id[a]) a++;
      while(b>0 && b<trg_pt_id.Dim() && trg_pt_id[b-1]==trg_pt_id[b]) b++;
      for(size_t i=a;i<b;i++){
        size_t pt_id=trg_pt_id[i]-trg_id_offset;
        trg_vel_final[pt_id*COORD_DIM+0]+=trg_vel[i*COORD_DIM+0];;
        trg_vel_final[pt_id*COORD_DIM+1]+=trg_vel[i*COORD_DIM+1];;
        trg_vel_final[pt_id*COORD_DIM+2]+=trg_vel[i*COORD_DIM+2];;
      }
    }
    trg_vel.Swap(trg_vel_final);
    pvfmm::Profile::Toc();
  }

  if(0){ // vis near-points
    pvfmm::Profile::Tic("Vis",&comm,true);
    typedef pvfmm::MPI_Node<Real_t> Node_t;
    typedef pvfmm::MPI_Tree<Node_t> Tree_t;
    typedef typename Node_t::NodeData NodeData_t;

    Tree_t pt_tree(comm);
    NodeData_t node_data;
    node_data.dim=COORD_DIM;
    node_data.max_depth=MAX_DEPTH;
    node_data.max_pts=10000000;

    { // Set node_data.pt_coord, node_data.pt_value
      Real_t scale_x, shift_x[COORD_DIM];
      {
        shift_x[0]=S_let.bbox[0];
        shift_x[1]=S_let.bbox[1];
        shift_x[2]=S_let.bbox[2];
        scale_x=S_let.bbox[3];
      }

      node_data.pt_coord=near_trg_coord;
      for(size_t i=0;i<node_data.pt_coord.Dim()/COORD_DIM;i++){
        node_data.pt_coord[i*COORD_DIM+0]=node_data.pt_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
        node_data.pt_coord[i*COORD_DIM+1]=node_data.pt_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
        node_data.pt_coord[i*COORD_DIM+2]=node_data.pt_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
      }

      const Vec_t& x=S_let.S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);

      std::vector<Real_t> value;
      size_t ves_id_offset;
      { // Scatter trg points to vesicles
        const Vec_t& x=S_let.S->getPosition();
        size_t N_ves = x.getNumSubs(); // Number of vesicles
        size_t M_ves = x.getStride(); // Points per vesicle
        assert(M_ves==x.getGridDim().first*x.getGridDim().second);
        { // Get ves_id_offset
          long long disp=0;
          long long size=N_ves;
          MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
          ves_id_offset=disp-size;
        }
      }
      for(size_t i=0;i<near_trg_cnt.Dim();i++){
        for(size_t j=0;j<near_trg_cnt[i];j++){
          value.push_back(ves_id_offset+i);
        }
      }
      node_data.pt_value=value;
    }

    pt_tree.Initialize(&node_data);
    pt_tree.Write2File("near_pts");
    pvfmm::Profile::Toc();
  }
  if(0){ // vis trg_vel
    pvfmm::Profile::Tic("Vis",&comm,true);
    typedef pvfmm::MPI_Node<Real_t> Node_t;
    typedef pvfmm::MPI_Tree<Node_t> Tree_t;
    typedef typename Node_t::NodeData NodeData_t;

    Tree_t pt_tree(comm);
    NodeData_t node_data;
    node_data.dim=COORD_DIM;
    node_data.max_depth=MAX_DEPTH;
    node_data.max_pts=10000000;

    { // Set node_data.pt_coord, node_data.pt_value
      Real_t scale_x, shift_x[COORD_DIM];
      {
        shift_x[0]=S_let.bbox[0];
        shift_x[1]=S_let.bbox[1];
        shift_x[2]=S_let.bbox[2];
        scale_x=S_let.bbox[3];
      }

      node_data.pt_coord=T;
      for(size_t i=0;i<node_data.pt_coord.Dim()/COORD_DIM;i++){
        node_data.pt_coord[i*COORD_DIM+0]=node_data.pt_coord[i*COORD_DIM+0]*scale_x+shift_x[0];
        node_data.pt_coord[i*COORD_DIM+1]=node_data.pt_coord[i*COORD_DIM+1]*scale_x+shift_x[1];
        node_data.pt_coord[i*COORD_DIM+2]=node_data.pt_coord[i*COORD_DIM+2]*scale_x+shift_x[2];
      }
      node_data.pt_value=trg_vel;
    }

    pt_tree.Initialize(&node_data);
    pt_tree.Write2File("trg_vel");
    pvfmm::Profile::Toc();
  }

  return &trg_vel[0];
}

template<typename Surf_t>
void NearSingular<Surf_t>::EvaluateNearSingular(Vec_t& T_vel){
  size_t omp_p=omp_get_max_threads();

  this->EvaluateNearSingular();

  const Vec_t& x=T_vel;
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  assert(M_ves==x.getGridDim().first*x.getGridDim().second);
  assert(N_ves*M_ves*COORD_DIM==trg_vel.Dim());

  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){ // Set tree pt data
    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;

    for(size_t i=a;i<b;i++){
      // read each component of x
      const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
      const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
      const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

      for(size_t j=0;j<M_ves;j++){
        xk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+0];
        yk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+1];
        zk[j]=trg_vel[(i*M_ves+j)*COORD_DIM+2];
      }
    }
  }
}


template <class Real_t>
struct QuadraticPatch{

  public:

  QuadraticPatch(){}

  QuadraticPatch(Real_t* x, int dof_){
    dof=dof_;
    coeff.resize(9*dof);
    for(size_t i=0;i<dof;i++){
      Real_t tmp[3*3];
      for(size_t j=0;j<3;j++){
        tmp[j*3+0]= x[(j*3+1)*dof+i];
        tmp[j*3+1]=(x[(j*3+2)*dof+i]-x[(j*3+0)*dof+i])*0.5;
        tmp[j*3+2]=(x[(j*3+2)*dof+i]+x[(j*3+0)*dof+i])*0.5-x[(j*3+1)*dof+i];
      }
      for(size_t j=0;j<3;j++){
        coeff[i*9+0*3+j]= tmp[1*3+j];
        coeff[i*9+1*3+j]=(tmp[2*3+j]-tmp[0*3+j])*0.5;
        coeff[i*9+2*3+j]=(tmp[2*3+j]+tmp[0*3+j])*0.5-tmp[1*3+j];
      }
    }
  }

  void eval(Real_t x, Real_t y, Real_t* val){
    Real_t x_[3]={1,x,x*x};
    Real_t y_[3]={1,y,y*y};
    for(size_t k=0;k<dof;k++) val[k]=0;
    for(size_t k=0;k<dof;k++){
      for(size_t i=0;i<3;i++){
        for(size_t j=0;j<3;j++){
          val[k]+=coeff[i+3*j+9*k]*x_[i]*y_[j];
        }
      }
    }
  }

  void grad(Real_t x, Real_t y, Real_t* val){
    Real_t x_[3]={1,x,x*x};
    Real_t y_[3]={1,y,y*y};
    Real_t x__[3]={0,1,2*x};
    Real_t y__[3]={0,1,2*y};
    for(size_t k=0;k<2*dof;k++) val[k]=0;
    for(size_t k=0;k<dof;k++){
      for(size_t i=0;i<3;i++){
        for(size_t j=0;j<3;j++){
          val[k+dof*0]+=coeff[i+3*j+9*k]*x__[i]*y_[j];
          val[k+dof*1]+=coeff[i+3*j+9*k]*x_[i]*y__[j];
        }
      }
    }
  }

  private:

  int dof;
  std::vector<Real_t> coeff;
};

template<typename Real_t>
static void LegPoly(const Real_t* x, size_t n, size_t q, Real_t* y){
  if(q>0) for(size_t i=0;i<n;i++) y[i]=1.0;
  if(q>1) for(size_t i=0;i<n;i++) y[n+i]=x[i];
  for(size_t j=2;j<q;j++){
    Real_t inv_j=1.0/j;
    for(size_t i=0;i<n;i++){
      y[j*n+i]=((2.0*j-1.0)*x[i]*y[(j-1)*n+i] - (j-1.0)*y[(j-2)*n+i])*inv_j;
    }
  }

  // Normalize
  for(size_t j=0;j<q;j++){
    for(size_t i=0;i<n;i++){
      y[j*n+i]*=sqrt(2.0*j+1.0);
    }
  }
}

template<typename Surf_t>
void NearSingular<Surf_t>::StokesNearSingular(Real_t r_near, const size_t* trg_cnt, const size_t* trg_dsp,  Real_t* trg_coord, Real_t* trg_veloc){
  const Vec_t &force=*force_single;

  typedef Vec_t Vec_t;
  typedef typename Vec_t::scalars_type Sca_t;
  size_t omp_p=omp_get_max_threads();

  const int force_dim=COORD_DIM;
  const int veloc_dim=COORD_DIM;

  const Vec_t& x=S->getPosition();
  size_t N_ves = x.getNumSubs(); // Number of vesicles
  size_t M_ves = x.getStride(); // Points per vesicle
  int k0_max=x.getGridDim().first;
  int k1_max=x.getGridDim().second;
  assert(M_ves==k0_max*k1_max);

  std::vector<Real_t> pole_quad;
  { // compute quadrature to find pole
    size_t p=x.getShOrder();

    //Gauss-Legendre quadrature nodes and weights
    std::vector<Real_t> x(p+1),w(p+1);
    cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &x[0], &w[0]);

    std::vector<Real_t> leg((p+1)*(p+1));
    LegPoly(&x[0],x.size(),p+1,&leg[0]);
    pole_quad.resize(p+1,0);
    for(size_t j=0;j<p+1;j++){
      for(size_t i=0;i<p+1;i++){
        pole_quad[i]+=leg[j*(p+1)+i]*sqrt(2.0*j+1.0);
      }
    }
    for(size_t i=0;i<p+1;i++){
      pole_quad[i]*=w[i]*0.25/p;
    }
  }

  size_t interp_deg=8;

  //#define INTERP_X(i) i                                    // Equi-spaced points
  #define INTERP_X(i) (i?1.0+(i-1.0)/(interp_deg-2):0.0)   // Clustered points
  assert(INTERP_X(0)==0);

  //#define INTERP_Y(x,j) pow(x,j)                   // Polynomial basis
  #define INTERP_Y(x,j) pow(1.0/(1.0+x),j+1) // Laurent polynomial

  pvfmm::Matrix<Real_t> M(interp_deg,interp_deg);
  { // matrix for computing interpolation coefficients
    for(size_t i=0;i<interp_deg;i++){
      Real_t x=INTERP_X(i);
      for(size_t j=0;j<interp_deg;j++){
        M[i][j]=INTERP_Y(x,j);
      }
    }
    M=M.pinv();
  }

  { // Compute trg_veloc. [[ At this point we can offload. ]]
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;

      PVFMMVec_t s_coord(M_ves*COORD_DIM);
      PVFMMVec_t s_force(M_ves*force_dim);
      PVFMMVec_t s_veloc(M_ves*veloc_dim);

      PVFMMVec_t pole_coord(2*COORD_DIM);
      PVFMMVec_t pole_veloc(2*veloc_dim);

      QuadraticPatch<Real_t> patch_coord;
      QuadraticPatch<Real_t> patch_veloc;
      PVFMMVec_t patch_coord_(3*3*COORD_DIM);
      PVFMMVec_t patch_veloc_(3*3*veloc_dim);

      for(size_t i=a;i<b;i++){ // loop over all vesicles
        { // Set s_coord, s_force, s_veloc
          // read each component of x
          const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

          // read each component of qforce
          const Real_t* fxk=force.getSubN_begin(i)+0*M_ves;
          const Real_t* fyk=force.getSubN_begin(i)+1*M_ves;
          const Real_t* fzk=force.getSubN_begin(i)+2*M_ves;
          assert(force_dim==3);

          // read each component of veloc
          const Real_t* vsxk=S_vel->getSubN_begin(i)+0*M_ves;
          const Real_t* vsyk=S_vel->getSubN_begin(i)+1*M_ves;
          const Real_t* vszk=S_vel->getSubN_begin(i)+2*M_ves;
          assert(veloc_dim==3);

          for(size_t j=0;j<M_ves;j++){
            s_coord[j*COORD_DIM+0]=xk[j];
            s_coord[j*COORD_DIM+1]=yk[j];
            s_coord[j*COORD_DIM+2]=zk[j];

            s_force[j*force_dim+0]=fxk[j];
            s_force[j*force_dim+1]=fyk[j];
            s_force[j*force_dim+2]=fzk[j];

            s_veloc[j*veloc_dim+0]=vsxk[j];
            s_veloc[j*veloc_dim+1]=vsyk[j];
            s_veloc[j*veloc_dim+2]=vszk[j];
          }
        }
        { // Set pole values: pole_coord, pole_veloc
          pole_coord.SetZero();
          pole_veloc.SetZero();
          for(size_t k0=0;k0<k0_max;k0++){
            for(size_t k1=0;k1<k1_max;k1++){
              size_t k=k1+k0*k1_max;
              pole_coord[0*COORD_DIM+0]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+0];
              pole_coord[0*COORD_DIM+1]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+1];
              pole_coord[0*COORD_DIM+2]+=pole_quad[k0_max-1-k0]*s_coord[k*COORD_DIM+2];
              pole_coord[1*COORD_DIM+0]+=pole_quad[         k0]*s_coord[k*COORD_DIM+0];
              pole_coord[1*COORD_DIM+1]+=pole_quad[         k0]*s_coord[k*COORD_DIM+1];
              pole_coord[1*COORD_DIM+2]+=pole_quad[         k0]*s_coord[k*COORD_DIM+2];

              pole_veloc[0*veloc_dim+0]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+0];
              pole_veloc[0*veloc_dim+1]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+1];
              pole_veloc[0*veloc_dim+2]+=pole_quad[k0_max-1-k0]*s_veloc[k*veloc_dim+2];
              pole_veloc[1*veloc_dim+0]+=pole_quad[         k0]*s_veloc[k*veloc_dim+0];
              pole_veloc[1*veloc_dim+1]+=pole_quad[         k0]*s_veloc[k*veloc_dim+1];
              pole_veloc[1*veloc_dim+2]+=pole_quad[         k0]*s_veloc[k*veloc_dim+2];
            }
          }
          //std::cout<<pole_coord<<'\n';
        }

        PVFMMVec_t t_coord(trg_cnt[i]*COORD_DIM, &trg_coord[trg_dsp[i]*COORD_DIM], false);
        PVFMMVec_t t_veloc(trg_cnt[i]*veloc_dim, &trg_veloc[trg_dsp[i]*veloc_dim], false);
        t_veloc.SetZero();

        { // Subtract wrong near potential
          stokes_sl(&s_coord[0], M_ves, &s_force[0], 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
          for(size_t j=0;j<t_veloc.Dim();j++) t_veloc[j]=-t_veloc[j];
        }

        // Add corrected near potential
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          Real_t  t_coord_j[COORD_DIM];
          Real_t* t_veloc_j;

          { // Set t_coord, t_veloc
            t_coord_j[0]=t_coord[j*COORD_DIM+0];
            t_coord_j[1]=t_coord[j*COORD_DIM+1];
            t_coord_j[2]=t_coord[j*COORD_DIM+2];
            t_veloc_j=&t_veloc[j*veloc_dim];
          }

          { // Set t_veloc_j for surface mesh points and {continue;}
            bool surface_point=false;
            for(size_t k0=0;k0<k0_max;k0++){
              for(size_t k1=0;k1<k1_max;k1++){
                size_t k=k1+k0*k1_max;
                // Add self interaction
                if(t_coord_j[0]==s_coord[k*COORD_DIM+0])
                if(t_coord_j[1]==s_coord[k*COORD_DIM+1])
                if(t_coord_j[2]==s_coord[k*COORD_DIM+2]){
                  t_veloc_j[0]+=s_veloc[k*veloc_dim+0];
                  t_veloc_j[1]+=s_veloc[k*veloc_dim+1];
                  t_veloc_j[2]+=s_veloc[k*veloc_dim+2];
                  k0=k0_max; k1=k1_max;
                  surface_point=true;
                  break;
                }
              }
            }
            if(surface_point) continue;
          }

          { // Find nearest point on mesh and create patch
            int k0_, k1_; // mesh coordinates for nearest point
            { // Find nearest point on mesh
              Real_t r2_min=-1.0;
              for(size_t k0=0;k0<k0_max;k0++){
                for(size_t k1=0;k1<k1_max;k1++){
                  size_t k=k1+k0*k1_max;
                  Real_t dx=(t_coord[j*COORD_DIM+0]-s_coord[k*COORD_DIM+0]);
                  Real_t dy=(t_coord[j*COORD_DIM+1]-s_coord[k*COORD_DIM+1]);
                  Real_t dz=(t_coord[j*COORD_DIM+2]-s_coord[k*COORD_DIM+2]);
                  Real_t r2=dx*dx+dy*dy+dz*dz;
                  if(r2<r2_min || r2_min<0){
                    r2_min=r2;
                    k0_=k0;
                    k1_=k1;
                  }
                }
              }
            }
            { // Create patch
              { // (0,0)
                if(k0_>0){
                  int k=((k1_+k1_max-1)%k1_max)+(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+0)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+0)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+0)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+0)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+0)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+0)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }
              { // (0,1)
                if(k0_>0){
                  int k=  k1_                  +(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+1)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+1)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+1)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+1)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+1)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+1)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }
              { // (0,2)
                if(k0_>0){
                  int k=((k1_       +1)%k1_max)+(k0_-1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(0*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(0*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(0*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(0*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(0*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(0*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(0*3+2)*COORD_DIM+0]=pole_coord[0*COORD_DIM+0];
                  patch_coord_[(0*3+2)*COORD_DIM+1]=pole_coord[0*COORD_DIM+1];
                  patch_coord_[(0*3+2)*COORD_DIM+2]=pole_coord[0*COORD_DIM+2];

                  patch_veloc_[(0*3+2)*veloc_dim+0]=pole_veloc[0*veloc_dim+0];
                  patch_veloc_[(0*3+2)*veloc_dim+1]=pole_veloc[0*veloc_dim+1];
                  patch_veloc_[(0*3+2)*veloc_dim+2]=pole_veloc[0*veloc_dim+2];
                }
              }

              { // (1,0)
                int k=((k1_+k1_max-1)%k1_max)+k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }
              { // (1,1)
                int k=  k1_                  +k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }
              { // (1,2)
                int k=((k1_       +1)%k1_max)+k0_*k1_max;
                assert(k>=0 && k<M_ves);
                patch_coord_[(1*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                patch_coord_[(1*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                patch_coord_[(1*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                patch_veloc_[(1*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                patch_veloc_[(1*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                patch_veloc_[(1*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
              }

              { // (2,0)
                if(k0_<k0_max-1){
                  int k=((k1_+k1_max-1)%k1_max)+(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+0)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+0)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+0)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+0)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+0)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+0)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+0)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+0)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+0)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+0)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+0)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+0)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }
              { // (2,1)
                if(k0_<k0_max-1){
                  int k=  k1_                  +(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+1)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+1)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+1)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+1)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+1)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+1)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+1)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+1)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+1)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+1)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+1)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+1)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }
              { // (2,2)
                if(k0_<k0_max-1){
                  int k=((k1_       +1)%k1_max)+(k0_+1)*k1_max;
                  assert(k>=0 && k<M_ves);
                  patch_coord_[(2*3+2)*COORD_DIM+0]=s_coord[k*COORD_DIM+0];
                  patch_coord_[(2*3+2)*COORD_DIM+1]=s_coord[k*COORD_DIM+1];
                  patch_coord_[(2*3+2)*COORD_DIM+2]=s_coord[k*COORD_DIM+2];

                  patch_veloc_[(2*3+2)*veloc_dim+0]=s_veloc[k*veloc_dim+0];
                  patch_veloc_[(2*3+2)*veloc_dim+1]=s_veloc[k*veloc_dim+1];
                  patch_veloc_[(2*3+2)*veloc_dim+2]=s_veloc[k*veloc_dim+2];
                }else{
                  patch_coord_[(2*3+2)*COORD_DIM+0]=pole_coord[1*COORD_DIM+0];
                  patch_coord_[(2*3+2)*COORD_DIM+1]=pole_coord[1*COORD_DIM+1];
                  patch_coord_[(2*3+2)*COORD_DIM+2]=pole_coord[1*COORD_DIM+2];

                  patch_veloc_[(2*3+2)*veloc_dim+0]=pole_veloc[1*veloc_dim+0];
                  patch_veloc_[(2*3+2)*veloc_dim+1]=pole_veloc[1*veloc_dim+1];
                  patch_veloc_[(2*3+2)*veloc_dim+2]=pole_veloc[1*veloc_dim+2];
                }
              }

              patch_coord=QuadraticPatch<Real_t>(&patch_coord_[0],COORD_DIM);
              patch_veloc=QuadraticPatch<Real_t>(&patch_veloc_[0],veloc_dim);
            }
          }

          std::vector<Real_t> interp_coord(interp_deg*COORD_DIM);
          std::vector<Real_t> interp_veloc(interp_deg*veloc_dim,0);
          { // Find nearest point on patch (first interpolation point)
            Real_t sc[COORD_DIM];
            sc[0]=patch_coord_[(3*1+1)*COORD_DIM+0];
            sc[1]=patch_coord_[(3*1+1)*COORD_DIM+1];
            sc[2]=patch_coord_[(3*1+1)*COORD_DIM+2];

            Real_t x=0,y=0;
            while(true){
              Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                                    t_coord_j[1]-sc[1],
                                    t_coord_j[2]-sc[2]};
              Real_t dR2=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];

              Real_t sgrad[2*COORD_DIM];
              patch_coord.grad(x,y,sgrad);

              Real_t dxdx=sgrad[0]*sgrad[0]+sgrad[1]*sgrad[1]+sgrad[2]*sgrad[2];
              Real_t dydy=sgrad[3]*sgrad[3]+sgrad[4]*sgrad[4]+sgrad[5]*sgrad[5];
              Real_t dxdy=sgrad[0]*sgrad[3]+sgrad[1]*sgrad[4]+sgrad[2]*sgrad[5];
              Real_t dxdR=sgrad[0]*   dR[0]+sgrad[1]*   dR[1]+sgrad[2]*   dR[2];
              Real_t dydR=sgrad[3]*   dR[0]+sgrad[4]*   dR[1]+sgrad[5]*   dR[2];
              if(dxdR*dxdR/dxdx+dydR*dydR/dydy < dR2*1e-6) break;

              Real_t dx, dy;
              if(dxdy){
                dx=(dxdR/dxdy-dydR/dydy)/(dxdx/dxdy-dxdy/dydy);
                dy=(dydR/dxdy-dxdR/dxdx)/(dydy/dxdy-dxdy/dxdx);
              }else{
                dx=dxdR/dxdx;
                dy=dydR/dydy;
              }

              while(true){
                patch_coord.eval(x+dx,y+dy,sc);
                Real_t dR_[COORD_DIM]={t_coord_j[0]-sc[0],
                                       t_coord_j[1]-sc[1],
                                       t_coord_j[2]-sc[2]};
                Real_t dR2_=dR_[0]*dR_[0]+dR_[1]*dR_[1]+dR_[2]*dR_[2];
                if(dR2_<dR2){
                  x+=dx*2.0/3.0; y+=dy*2.0/3.0;
                  patch_coord.eval(x,y,sc);
                  break;
                }else{
                  dx*=0.5;
                  dy*=0.5;
                }
              }
            }

            // Set first interpolation point
            interp_coord[0*COORD_DIM+0]=sc[0];
            interp_coord[0*COORD_DIM+1]=sc[1];
            interp_coord[0*COORD_DIM+2]=sc[2];
            patch_veloc.eval(x,y,&interp_veloc[0*veloc_dim]);

            // For visualization of near point projections
            //t_coord[j*COORD_DIM+0]=sc[0];
            //t_coord[j*COORD_DIM+1]=sc[1];
            //t_coord[j*COORD_DIM+2]=sc[2];
          }

          { // Interpolate and find target velocity t_veloc_j
            Real_t dR[COORD_DIM]={t_coord_j[0]-interp_coord[0],
                                t_coord_j[1]-interp_coord[1],
                                t_coord_j[2]-interp_coord[2]};
            Real_t dR_norm=sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
            Real_t OOdR=1.0/dR_norm;

            for(size_t l=1;l<interp_deg;l++){
              Real_t x=INTERP_X(l);
              interp_coord[l*COORD_DIM+0]=interp_coord[0*COORD_DIM+0]+dR[0]*OOdR*r_near*x;
              interp_coord[l*COORD_DIM+1]=interp_coord[0*COORD_DIM+1]+dR[1]*OOdR*r_near*x;
              interp_coord[l*COORD_DIM+2]=interp_coord[0*COORD_DIM+2]+dR[2]*OOdR*r_near*x;
            }

            // Compute velocity at interpolation points
            stokes_sl(&s_coord[0], M_ves, &s_force[0], 1, &interp_coord[1*COORD_DIM], interp_deg-1, &interp_veloc[1*veloc_dim], NULL);

            // Interpolate
            for(size_t k=0;k<veloc_dim;k++){
              pvfmm::Matrix<Real_t> y(interp_deg,1);
              pvfmm::Matrix<Real_t> x(interp_deg,1);
              pvfmm::Matrix<Real_t> coeff(interp_deg,1);

              Real_t x_=dR_norm/r_near;
              for(size_t l=0;l<interp_deg;l++){
                y[l][0]=interp_veloc[l*veloc_dim+k];
                x[l][0]=INTERP_Y(x_,l);
              }

              coeff=M*y;
              for(size_t l=0;l<interp_deg;l++) t_veloc_j[k]+=coeff[l][0]*x[l][0];
            }
          }
        }
      }
    }
  }
}



template<typename Real_t>
static void u_ref(const Real_t* coord, int n, Real_t* out){ //Analytical velocity for sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    Real_t r_2=0;
    r_2+=c[0]*c[0];
    r_2+=c[1]*c[1];
    r_2+=c[2]*c[2];
    Real_t r=sqrt(r_2);

    Real_t cos_t=c[0]/r;
    Real_t sin_t=sqrt(1-cos_t*cos_t);
    Real_t R0_r=R0/r;
    Real_t ur= cos_t*(1.0-1.50*R0_r+0.50*R0_r*R0_r*R0_r);
    Real_t ut=-sin_t*(1.0-0.75*R0_r-0.25*R0_r*R0_r*R0_r);
    out[i*COORD_DIM+0]=cos_t*ur-sin_t*ut-1.0;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*COORD_DIM+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*COORD_DIM+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;
  }
}

template<typename Real_t>
static void force(const Real_t* coord, int n, Real_t* out){ // Force on sphere
  Real_t R0=1.0;
  for(int i=0;i<n;i++){
    const Real_t* c=&coord[i*COORD_DIM];
    Real_t r_2=0;
    r_2+=c[0]*c[0];
    r_2+=c[1]*c[1];
    r_2+=c[2]*c[2];
    Real_t r=sqrt(r_2);

    Real_t cos_t=c[0]/r;
    Real_t sin_t=sqrt(1-cos_t*cos_t);
    Real_t R0_r=R0/r;
    Real_t ur=-1.5*cos_t;
    Real_t ut=+1.5*sin_t;
    out[i*COORD_DIM+0]=cos_t*ur-sin_t*ut;

    Real_t r_yz=sqrt(c[1]*c[1]+c[2]*c[2]);
    out[i*COORD_DIM+1]=(sin_t*ur+cos_t*ut)*c[1]/r_yz;
    out[i*COORD_DIM+2]=(sin_t*ur+cos_t*ut)*c[2]/r_yz;
  }
}

template<typename Surf_t>
void NearSingular<Surf_t>::TestNearSingular(){
  typedef typename Vec_t::scalars_type Sca_t;
  typedef typename Vec_t::array_type Arr_t;
  typedef OperatorsMats<Arr_t> Mats_t;
  typedef VesInteraction<Real_t> Interaction_t;
  typedef InterfacialVelocity<Surf_t, Interaction_t> IntVel_t;

  int rank, size;
  MPI_Comm comm=MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  size_t nVec=2/size;
  if(!nVec) nVec=1;
  srand48(rank);

  //Set parameters
  Parameters<Real_t> sim_par;
  sim_par.sh_order = 32;
  sim_par.rep_up_freq = 32;

  //Create vectors
  Vec_t x0(nVec, sim_par.sh_order); // coordinates
  Vec_t v0(nVec, sim_par.sh_order); // velocity
  Vec_t v0_ref(nVec, sim_par.sh_order); // reference velocity
  Vec_t f0(nVec, sim_par.sh_order); // force
  int fLen = x0.getStride();
  { //Set coordinate values
    int imax(x0.getGridDim().first);
    int jmax(x0.getGridDim().second);
    Real_t scal=pow((Real_t)nVec*size,1.0/3.0);

    std::vector<Real_t> qx;
    { // compute legendre node points
      size_t p=x0.getShOrder();
      qx.resize(p+1);
      std::vector<Real_t> qw(p+1);
      cgqf(p+1, 1, 0.0, 0.0, -1.0, 1.0, &qx[0], &qw[0]);
    }

    for(size_t k=0;k<nVec;k++){
      Real_t* x_k=x0.getSubN_begin(k)+0*fLen;
      Real_t* y_k=x0.getSubN_begin(k)+1*fLen;
      Real_t* z_k=x0.getSubN_begin(k)+2*fLen;

      Real_t* vx_k=v0.getSubN_begin(k)+0*fLen;
      Real_t* vy_k=v0.getSubN_begin(k)+1*fLen;
      Real_t* vz_k=v0.getSubN_begin(k)+2*fLen;
      Real_t* fx_k=f0.getSubN_begin(k)+0*fLen;
      Real_t* fy_k=f0.getSubN_begin(k)+1*fLen;
      Real_t* fz_k=f0.getSubN_begin(k)+2*fLen;

      Real_t coord[3];
      coord[0]=scal*drand48();
      coord[1]=scal*drand48();
      coord[2]=scal*drand48();
      for(size_t i=0;i<imax;i++){
        Real_t cos_t=qx[i];
        //Real_t cos_t=cos((i+1)*M_PI/(imax+1));
        Real_t sin_t=sqrt(1.0-cos_t*cos_t);
        for(size_t j=0;j<jmax;j++){

          if(0){
            x_k[j+i*jmax]=cos_t;
            y_k[j+i*jmax]=sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=sin_t*cos(j*2*M_PI/jmax);
            if(k){
              x_k[j+i*jmax]*=0.49;
              y_k[j+i*jmax]*=0.49;
              z_k[j+i*jmax]*=0.49;
              y_k[j+i*jmax]+=1.49+1e-4;
            }
          }else{
            x_k[j+i*jmax]=coord[0]+scal*cos_t;
            y_k[j+i*jmax]=coord[1]+scal*sin_t*sin(j*2*M_PI/jmax);
            z_k[j+i*jmax]=coord[2]+scal*sin_t*cos(j*2*M_PI/jmax);
          }

          vx_k[j+i*jmax]=0;
          vy_k[j+i*jmax]=0;
          vz_k[j+i*jmax]=0;
          fx_k[j+i*jmax]=0;
          fy_k[j+i*jmax]=0;
          fz_k[j+i*jmax]=0;

          if(!rank && k==0){
            Real_t c[3]={x_k[j+i*jmax],
                       y_k[j+i*jmax],
                       z_k[j+i*jmax]};
            Real_t f[3]={0,0,0};
            force(c,1,f);
            fx_k[j+i*jmax]=f[0];
            fy_k[j+i*jmax]=f[1];
            fz_k[j+i*jmax]=f[2];
          }
        }
      }
    }
  }

  //Reading operators from file
  bool readFromFile = true;
  Mats_t mats(readFromFile, sim_par);

  //Creating objects
  COUT("Creating the surface object");
  Surf_t S(x0, mats);

  //Setting the background flow
  ShearFlow<Vec_t> vInf(sim_par.bg_flow_param);

  Interaction_t interaction(NULL);
  IntVel_t F(S, interaction, mats, sim_par, vInf);

  // Compute self interaction
  Vec_t veloc_self(nVec, sim_par.sh_order);
  F.stokes(f0,veloc_self);

  Vec_t qforce(nVec,sim_par.sh_order);
  { //Incorporating the quadrature weights and area into the force
    Sca_t quad_weights_;
    { // quadrature weights
      quad_weights_.resize(1,sim_par.sh_order);
      quad_weights_.getDevice().Memcpy(quad_weights_.begin(),
          mats.quad_weights_,
          quad_weights_.size() * sizeof(Real_t),
          Surf_t::device_type::MemcpyDeviceToDevice);
    }
    xv(S.getAreaElement(), f0, qforce);
    ax<Sca_t>(quad_weights_, qforce, qforce);
  }

  NearSingular<Surf_t> near_singular(comm);
  near_singular.SetSrcCoord(S);
  near_singular.SetTrgCoord(S);
  near_singular.SetSurfaceVel(veloc_self);
  near_singular.SetDensity(&qforce, NULL);

  pvfmm::Profile::Tic("NearInteractions",&comm);
  near_singular.EvaluateNearSingular();
  pvfmm::Profile::Toc();
}

#endif
