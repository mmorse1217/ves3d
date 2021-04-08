#include <omp.h>
#include <iostream>

#include <ompUtils.h>
#include <parUtils.h>
#include <profile.hpp>
#include <mortonid.hpp>
#include <SphericalHarmonics.h>
#include <mpi_tree.hpp> // Only for vis

#define VES_STRIDE (2+(2*sh_order_)*(1+sh_order_))

template<typename Real_t>
NearSingular<Real_t>::NearSingular(Real_t box_size, Real_t repul_dist, MPI_Comm c){
  box_size_=box_size;
  repul_dist_=repul_dist;
  comm=c;
  S=NULL;
  qforce_single=NULL;
  qforce_double=NULL;
  force_double=NULL;
  S_vel=NULL;

  update_setup =update_setup  | NearSingular::UpdateSrcCoord;
  update_setup =update_setup  | NearSingular::UpdateTrgCoord;

  update_direct=update_direct | NearSingular::UpdateSrcCoord;
  update_direct=update_direct | NearSingular::UpdateDensitySL;
  update_direct=update_direct | NearSingular::UpdateDensityDL;
  update_direct=update_direct | NearSingular::UpdateTrgCoord;

  update_interp=update_interp | NearSingular::UpdateSrcCoord;
  update_interp=update_interp | NearSingular::UpdateSurfaceVel;
  update_interp=update_interp | NearSingular::UpdateDensitySL;
  update_interp=update_interp | NearSingular::UpdateDensityDL;
  update_interp=update_interp | NearSingular::UpdateTrgCoord;
}


template<typename Real_t>
void NearSingular<Real_t>::SetSrcCoord(const PVFMMVec_t& src_coord, int sh_order){
  update_direct=update_direct | NearSingular::UpdateSrcCoord;
  update_interp=update_interp | NearSingular::UpdateSrcCoord;
  update_setup =update_setup  | NearSingular::UpdateSrcCoord;
  sh_order_=sh_order;
  S=&src_coord;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetSurfaceVel(const PVFMMVec_t* S_vel_ptr){
  update_interp=update_interp | NearSingular::UpdateSurfaceVel;
  S_vel=S_vel_ptr;
}

template<typename Real_t>
void NearSingular<Real_t>::SetDensitySL(const PVFMMVec_t* qforce_single_){
  update_direct=update_direct | NearSingular::UpdateDensitySL;
  update_interp=update_interp | NearSingular::UpdateDensitySL;
  qforce_single=qforce_single_;
}

template<typename Real_t>
void NearSingular<Real_t>::SetDensityDL(const PVFMMVec_t* qforce_double_, const PVFMMVec_t* force_double_){
  update_direct=update_direct | NearSingular::UpdateDensityDL;
  update_interp=update_interp | NearSingular::UpdateDensityDL;
  qforce_double=qforce_double_;
  force_double=force_double_;
}


template<typename Real_t>
void NearSingular<Real_t>::SetTrgCoord(Real_t* trg_coord, size_t N, bool trg_is_surf_){ // TODO: change to const Real_t*
  update_direct=update_direct | NearSingular::UpdateTrgCoord;
  update_interp=update_interp | NearSingular::UpdateTrgCoord;
  update_setup =update_setup  | NearSingular::UpdateTrgCoord;
  trg_is_surf=trg_is_surf_;
  T.ReInit(N*COORD_DIM,trg_coord);
}

template<typename Real_t>
void NearSingular<Real_t>::SetupCoordData(){
  INFO("near coord setup");
  assert(S);
  Real_t near=2.0/sqrt((Real_t)sh_order_); // TODO: some function of sh_order and accuracy
  if(!(update_setup & (NearSingular::UpdateSrcCoord | NearSingular::UpdateTrgCoord))) return;
  update_setup=update_setup & ~(NearSingular::UpdateSrcCoord | NearSingular::UpdateTrgCoord);

  int np, rank;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&rank);
  size_t omp_p=omp_get_max_threads();

  pvfmm::Profile::Tic("NearSetup",&comm,true);
  bool prof_state=pvfmm::Profile::Enable(false);
  struct{
    pvfmm::Vector<pvfmm::MortonId> mid; // MortonId of leaf nodes
    pvfmm::Vector<size_t> pt_cnt;       // Point count
    pvfmm::Vector<size_t> pt_dsp;       // Point displ
    pvfmm::Vector<pvfmm::MortonId> mins;// First non-ghost node

    PVFMMVec_t pt_coord;     // All point coordinates
    pvfmm::Vector<size_t> pt_vesid;     // = pt_id/M_ves
    pvfmm::Vector<size_t> pt_id;        // Scatter id
  } S_let;
  { // Construct S_let
    pvfmm::Profile::Tic("VesLET",&comm,true);

    size_t M_ves = VES_STRIDE;                 // Points per vesicle
    size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
    assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

    size_t tree_depth=0;
    { // Determine bbox, r_near, tree_depth
      pvfmm::Profile::Tic("PtData",&comm,true);
      Real_t* bbox=coord_setup.bbox;
      Real_t& r_near=coord_setup.r_near;

      Real_t r_ves=0;
      { // Determine r_ves
        std::vector<Real_t> r2_ves_(omp_p);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=((tid+0)*N_ves)/omp_p;
          size_t b=((tid+1)*N_ves)/omp_p;

          Real_t r2_ves=0;
          Real_t one_over_M=1.0/M_ves;
          for(size_t i=a;i<b;i++){ // compute r2_ves
            const Real_t* Si=&S[0][i*M_ves*COORD_DIM];
            Real_t center_coord[COORD_DIM]={0,0,0};
            for(size_t j=0;j<M_ves;j++){
              center_coord[0]+=Si[j*COORD_DIM+0];
              center_coord[1]+=Si[j*COORD_DIM+1];
              center_coord[2]+=Si[j*COORD_DIM+2];
            }
            center_coord[0]*=one_over_M;
            center_coord[1]*=one_over_M;
            center_coord[2]*=one_over_M;
            for(size_t j=0;j<M_ves;j++){
              Real_t dx=(Si[j*COORD_DIM+0]-center_coord[0]);
              Real_t dy=(Si[j*COORD_DIM+1]-center_coord[1]);
              Real_t dz=(Si[j*COORD_DIM+2]-center_coord[2]);
              Real_t r2=dx*dx+dy*dy+dz*dz;
              r2_ves=std::max(r2_ves,r2);
            }
          }
          r2_ves_[tid]=r2_ves;
        }

        { // Determine r_ves (global max)
          double r_ves_loc=0, r_ves_glb=0;
          for(size_t tid=0;tid<omp_p;tid++){
            r_ves_loc=std::max(r2_ves_[tid], r_ves_loc);
          }
          r_ves_loc=sqrt(r_ves_loc);
          MPI_Allreduce(&r_ves_loc, &r_ves_glb, 1, MPI_DOUBLE, MPI_MAX, comm);
          r_ves=r_ves_glb;
        }
      }

      r_near=r_ves*near; // r_near is some function of r_ves.
      if(box_size_>0 && 2*r_near+r_ves>box_size_){ // domain too small; abort
        COUTDEBUG("Domain too small for vesicle size. Multiple copies of a point can be NEAR a vesicle.");
        assert(false);
        exit(0);
      }

      if(box_size_<=0){ // Determine bbox, tree_depth
        Real_t scale_x, shift_x[COORD_DIM];
        Real_t scale_tmp;
        { // determine bounding box
          Real_t s0, x0[COORD_DIM];
          Real_t s1, x1[COORD_DIM];
          PVFMMBoundingBox(      N_ves*M_ves, &S[0][0], &s0, x0, comm);
          PVFMMBoundingBox(T.Dim()/COORD_DIM, &   T[0], &s1, x1, comm);

          Real_t c0[COORD_DIM]={(0.5-x0[0])/s0, (0.5-x0[1])/s0, (0.5-x0[2])/s0};
          Real_t c1[COORD_DIM]={(0.5-x1[0])/s1, (0.5-x1[1])/s1, (0.5-x1[2])/s1};

          scale_tmp=0;
          scale_tmp=std::max(scale_tmp, fabs(c0[0]-c1[0]));
          scale_tmp=std::max(scale_tmp, fabs(c0[1]-c1[1]));
          scale_tmp=std::max(scale_tmp, fabs(c0[2]-c1[2]));
          scale_tmp=1.0/(scale_tmp+1/s0+1/s1);

          shift_x[0]=0.5-(c0[0]+c1[0])*scale_tmp/2.0;
          shift_x[1]=0.5-(c0[1]+c1[1])*scale_tmp/2.0;
          shift_x[2]=0.5-(c0[2]+c1[2])*scale_tmp/2.0;
        }

        { // scale_x, pt_tree_depth
          if(scale_tmp==0){
            scale_tmp=1.0;
            r_near=1.0;
          }
          Real_t domain_length=1.0/scale_tmp+4*r_near;
          Real_t leaf_length=r_near;
          scale_x=1.0/leaf_length;
          while(domain_length*scale_x>1.0 && tree_depth<MAX_DEPTH-1){
            scale_x*=0.5;
            tree_depth++;
          }
        }
        for(size_t j=0;j<COORD_DIM;j++){ // Update shift_x
          shift_x[j]=((shift_x[j]/scale_tmp)+2*r_near)*scale_x;
        }
        coord_setup.bbox[0]=shift_x[0];
        coord_setup.bbox[1]=shift_x[1];
        coord_setup.bbox[2]=shift_x[2];
        coord_setup.bbox[3]=scale_x;
      }else{
        coord_setup.bbox[0]=0;
        coord_setup.bbox[1]=0;
        coord_setup.bbox[2]=0;
        coord_setup.bbox[3]=1.0/box_size_;

        // Determine tree depth
        Real_t leaf_size=1.0/coord_setup.bbox[3];
        while(leaf_size*0.5>r_near && tree_depth<MAX_DEPTH-1){
          leaf_size*=0.5;
          tree_depth++;
        }
      }
      pvfmm::Profile::Toc();
    }

    { // Construct local tree S_let
      pvfmm::Profile::Tic("LocalTree",&comm,true);
      PVFMMVec_t           & pt_coord=S_let.pt_coord;
      pvfmm::Vector<size_t>& pt_vesid=S_let.pt_vesid;
      pvfmm::Vector<size_t>& pt_id   =S_let.pt_id   ;

      pvfmm::Vector<pvfmm::MortonId>& let_mid   =S_let.mid;
      pvfmm::Vector<size_t>&          let_pt_cnt=S_let.pt_cnt;
      pvfmm::Vector<size_t>&          let_pt_dsp=S_let.pt_dsp;
      pvfmm::Vector<pvfmm::MortonId>& let_mins  =S_let.mins;

      { // build scatter-indices (pt_id) and tree (let_mid, let_pt_cnt, let_pt_dsp)
        pvfmm::Vector<pvfmm::MortonId> pt_mid(N_ves*M_ves);
        { // build pt_mid
          Real_t scale_x, shift_x[COORD_DIM];
          { // set scale_x, shift_x
            shift_x[0]=coord_setup.bbox[0];
            shift_x[1]=coord_setup.bbox[1];
            shift_x[2]=coord_setup.bbox[2];
            scale_x=coord_setup.bbox[3];
          }

          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){
            Real_t c[COORD_DIM];
            size_t a=((tid+0)*N_ves*M_ves)/omp_p;
            size_t b=((tid+1)*N_ves*M_ves)/omp_p;
            for(size_t i=a;i<b;i++){
              for(size_t k=0;k<COORD_DIM;k++){
                c[k]=S[0][i*COORD_DIM+k]*scale_x+shift_x[k];
                while(c[k]< 0.0) c[k]+=1.0;
                while(c[k]>=1.0) c[k]-=1.0;
              }
              pt_mid[i]=pvfmm::MortonId(c,tree_depth);
            }
          }
        }

        pt_id   .ReInit(N_ves*M_ves);
        pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm);
        pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
        { // build let_mins
          let_mins.ReInit(np);
          MPI_Allgather(&  pt_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                        &let_mins[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
          if(rank) assert(let_mins[rank]!=let_mins[rank-1]);
          let_mins[0]=pvfmm::MortonId(0,0,0,tree_depth);
        }
        { // Exchange shared octant with neighbour
          int send_size=0;
          int recv_size=0;
          if(rank<np-1){ // send_size
            send_size=pt_mid.Dim()-(std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mins[rank+1])-&pt_mid[0]);
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
            for(size_t i=0;i<recv_size;i++) pt_mid_new[i]=let_mins[rank];
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

        { // set let_mid
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
        }
        { // set let_pt_dsp
          let_pt_dsp.ReInit(let_mid.Dim());
          #pragma omp parallel for
          for(size_t i=0;i<let_mid.Dim();i++){
            let_pt_dsp[i]=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mid[i])-&pt_mid[0];
          }
        }
        { // set let_pt_cnt
          let_pt_cnt.ReInit(let_mid.Dim());
          #pragma omp parallel for
          for(size_t i=1;i<let_mid.Dim();i++){
            let_pt_cnt[i-1]=let_pt_dsp[i]-let_pt_dsp[i-1];
          }
          if(let_mid.Dim()) let_pt_cnt[let_mid.Dim()-1]=pt_mid.Dim()-let_pt_dsp[let_mid.Dim()-1];
        }
      }
      { // scatter pt_coord
        pt_coord=S[0];
        pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
      }
      { // scatter pt_vesid
        size_t ves_id_offset;
        { // Get ves_id_offset
          long long disp=0;
          long long size=N_ves;
          MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
          ves_id_offset=disp-size;
        }

        pt_vesid.ReInit(N_ves*M_ves);
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){ // set pt_vesid
          size_t a=((tid+0)*N_ves)/omp_p;
          size_t b=((tid+1)*N_ves)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t j=0;j<M_ves;j++){
              pt_vesid[i*M_ves+j]=ves_id_offset+i;
            }
          }
        }
        pvfmm::par::ScatterForward(pt_vesid, pt_id, comm);
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
            pvfmm::Vector<pvfmm::MortonId>&  let_mid=S_let. mid;
            pvfmm::Vector<pvfmm::MortonId>& let_mins=S_let.mins;

            Real_t coord[COORD_DIM];
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
                if(box_size_>0){
                  for(size_t k=0;k<COORD_DIM;k++){
                    if(c[k]< 0.0) c[k]+=1.0;
                    if(c[k]>=1.0) c[k]-=1.0;
                  }
                }
                if(c[0]>=0 && c[1]>=0 && c[2]>=0 && c[0]<1.0 && c[1]<1.0 && c[2]<1.0){
                  pvfmm::MortonId m(c,tree_depth+1);
                  if((rank && m<let_mins[rank]) || (rank<np-1 && m>=let_mins[rank+1])){
                    int pid=std::lower_bound(&let_mins[0], &let_mins[0]+np, m)-&let_mins[0]-1;
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
        { // Set snode_cnt, snode_id
          snode_cnt.SetZero();
          snode_id.ReInit(shared_pair.Dim());
          // TODO: parallelize this loop
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

        { // Sort snode_id for each process
          #pragma omp parallel for
          for(size_t tid=0;tid<omp_p;tid++){
            size_t a=(tid+0)*shared_pair.Dim()/omp_p;
            size_t b=(tid+1)*shared_pair.Dim()/omp_p;
            size_t pid_a=(a<shared_pair.Dim()?shared_pair[a].key:np);
            size_t pid_b=(b<shared_pair.Dim()?shared_pair[b].key:np);
            for(size_t i=pid_a;i<pid_b;i++){
              a=(i+0<np?snode_dsp[i+0]:shared_pair.Dim());
              b=(i+1<np?snode_dsp[i+1]:shared_pair.Dim());
              if(a<b) std::sort(&snode_id[0]+a,&snode_id[0]+b);
            }
          }
        }
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

  { // find vesicle, near target pairs TODO: cleanup
    pvfmm::Profile::Tic("TrgNear",&comm,true);

    PVFMMVec_t&            near_trg_coord  =coord_setup.near_trg_coord  ;
    pvfmm::Vector<size_t>& near_trg_cnt    =coord_setup.near_trg_cnt    ;
    pvfmm::Vector<size_t>& near_trg_dsp    =coord_setup.near_trg_dsp    ;
    pvfmm::Vector<size_t>& near_trg_scatter=coord_setup.near_trg_scatter;
    pvfmm::Vector<size_t>& near_trg_pt_id  =coord_setup.near_trg_pt_id  ;
    pvfmm::Vector<size_t>& near_ves_pt_id  =coord_setup.near_ves_pt_id  ;

    pvfmm::Vector<long long> pt_vesid;
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
          shift_x[0]=coord_setup.bbox[0];
          shift_x[1]=coord_setup.bbox[1];
          shift_x[2]=coord_setup.bbox[2];
          scale_x=coord_setup.bbox[3];
        }
        assert(S_let.mid.Dim());
        size_t tree_depth=S_let.mid[0].GetDepth();
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          Real_t c[COORD_DIM];
          size_t a=((tid+0)*T.Dim()/COORD_DIM)/omp_p;
          size_t b=((tid+1)*T.Dim()/COORD_DIM)/omp_p;
          for(size_t i=a;i<b;i++){
            for(size_t k=0;k<COORD_DIM;k++){
              c[k]=pt_coord[i*COORD_DIM+k]*scale_x+shift_x[k];
              while(c[k]< 0.0) c[k]+=1.0;
              while(c[k]>=1.0) c[k]-=1.0;
            }
            pt_mid[i]=pvfmm::MortonId(c,tree_depth);
          }
        }
      }

      pt_vesid.ReInit(pt_coord.Dim()/COORD_DIM);
      if(trg_is_surf){ // Set pt_vesid
        size_t N_ves = S->Dim()/(VES_STRIDE*COORD_DIM); // Number of vesicles
        size_t M_ves = pt_vesid.Dim()/N_ves;
        size_t ves_id_offset;
        { // Get ves_id_offset
          long long disp=0;
          long long size=N_ves;
          MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
          ves_id_offset=disp-size;
        }
        assert(N_ves*M_ves==pt_vesid.Dim());
        for(size_t i=0;i<N_ves;i++)
        for(size_t j=0;j<M_ves;j++){
          pt_vesid[i*M_ves+j]=ves_id_offset+i;
        }
      }else{
        for(size_t i=0;i<pt_vesid.Dim();i++) pt_vesid[i]=-1;
      }

      { // Sort pt data (pt_mid, pt_coord, pt_id)
        pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm, &S_let.mins[rank]);
        pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
        pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
        pvfmm::par::ScatterForward(pt_vesid, pt_id, comm);
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

    { // Find near trg, ves points
      pvfmm::Profile::Tic("NearPair",&comm,true);
      std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_(omp_p); // (loc_trg_id, S_let.pt_id)
      #pragma omp parallel num_threads(omp_p)
      {
        size_t tid=omp_get_thread_num();
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_[tid];
        size_t tree_depth; Real_t r2_near;
        { // Set tree_depth, r_near
          tree_depth=S_let.mid[0].GetDepth();
          r2_near=coord_setup.r_near;
          r2_near*=r2_near;
        }
        Real_t s=pow(0.5,tree_depth);

        size_t FLOP=0;
        size_t a=((tid+0)*tree_mid.Dim())/omp_p;
        size_t b=((tid+1)*tree_mid.Dim())/omp_p;
        for(size_t i=a;i<b;i++){
          size_t tcnt=tree_pt_cnt[i];
          size_t tdsp=tree_pt_dsp[i];
          PVFMMVec_t tcoord;
          pvfmm::Vector<long long> tvesid;
          if(tcnt){ // Set t_coord
            tcoord.ReInit(tcnt*COORD_DIM,&pt_coord[tdsp*COORD_DIM],false);
            tvesid.ReInit(tcnt          ,&pt_vesid[tdsp          ],false);
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
              if(box_size_>0){
                for(size_t k=0;k<COORD_DIM;k++){
                  if(c[k]< 0.0) c[k]+=1.0;
                  if(c[k]>=1.0) c[k]-=1.0;
                }
              }
              if(c[0]>=0 && c[1]>=0 && c[2]>=0 && c[0]<1.0 && c[1]<1.0 && c[2]<1.0){
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
                Real_t dx=fabs(scoord[k][s*COORD_DIM+0]-tcoord[t*COORD_DIM+0]);
                Real_t dy=fabs(scoord[k][s*COORD_DIM+1]-tcoord[t*COORD_DIM+1]);
                Real_t dz=fabs(scoord[k][s*COORD_DIM+2]-tcoord[t*COORD_DIM+2]);
                if(box_size_>0){
                  while(dx>box_size_*0.5) dx-=box_size_;
                  while(dy>box_size_*0.5) dy-=box_size_;
                  while(dz>box_size_*0.5) dz-=box_size_;
                }
                Real_t r2=dx*dx+dy*dy+dz*dz;
                if(r2<r2_near && svesid[k][s]!=tvesid[t]){
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
      size_t M_ves = VES_STRIDE;                 // Points per vesicle
      size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
      assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

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

      { // Compute near_trg_cnt, near_trg_dsp
        near_trg_cnt.ReInit(N_ves);near_trg_cnt.SetZero();
        near_trg_dsp.ReInit(N_ves);near_trg_dsp[0]=0;

        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(near_trg_pt_id.Dim()*(tid+0))/omp_p;
          size_t b=(near_trg_pt_id.Dim()*(tid+1))/omp_p;
          for(size_t i=a;i<b;i++) near_ves_pt_id[i]-=ves_pt_id_offset;
        }
        #pragma omp parallel for
        for(size_t tid=0;tid<omp_p;tid++){
          size_t a=(near_trg_pt_id.Dim()*(tid+0))/omp_p;
          size_t b=(near_trg_pt_id.Dim()*(tid+1))/omp_p;
          while(a>0 && a<near_trg_pt_id.Dim() && near_ves_pt_id[a-1]/M_ves==near_ves_pt_id[a]/M_ves) a++;
          while(b>0 && b<near_trg_pt_id.Dim() && near_ves_pt_id[b-1]/M_ves==near_ves_pt_id[b]/M_ves) b++;
          for(size_t i=a;i<b;i++) near_trg_cnt[near_ves_pt_id[i]/M_ves]++;
        }
        pvfmm::omp_par::scan(&near_trg_cnt[0], &near_trg_dsp[0], N_ves);
        assert(near_trg_dsp[N_ves-1]+near_trg_cnt[N_ves-1]==near_trg_coord.Dim()/COORD_DIM);
      }
      pvfmm::Profile::Toc();
    }

    { // Determine near_trg_scatter, near_trg_pt_id
      pvfmm::par::SortScatterIndex(near_trg_pt_id, near_trg_scatter, comm, &trg_id_offset);
      pvfmm::par::ScatterForward  (near_trg_pt_id, near_trg_scatter, comm);
    }

    if(box_size_>0){ // periodic translation for target points
      PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
      pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
      pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;

      size_t M_ves = VES_STRIDE;                 // Points per vesicle
      size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
      assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

      for(size_t i=0;i<N_ves;i++){ // loop over all vesicles
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over near tagets
          size_t trg_idx=trg_dsp[i]+j;
          size_t src_idx=coord_setup.near_ves_pt_id[trg_idx];
          Real_t ves_c[COORD_DIM]={S[0][src_idx*COORD_DIM+0],
                                   S[0][src_idx*COORD_DIM+1],
                                   S[0][src_idx*COORD_DIM+2]};
          for(size_t k=0;k<COORD_DIM;k++){
            Real_t c=trg_coord[trg_idx*COORD_DIM+k];
            while(c-ves_c[k]> box_size_*0.5) c-=box_size_;
            while(c-ves_c[k]<-box_size_*0.5) c+=box_size_;
            trg_coord[trg_idx*COORD_DIM+k]=c;
          }
        }
      }
    }

    pvfmm::Profile::Toc();
  }
  { // projection
    pvfmm::Profile::Tic("Proj",&comm,true);
    PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
    pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
    pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;

    pvfmm::Vector<char>& is_extr_pt=coord_setup.is_extr_pt;
    PVFMMVec_t& proj_patch_param=coord_setup.proj_patch_param;
    PVFMMVec_t& proj_coord      =coord_setup.proj_coord      ;
    PVFMMVec_t& repl_force      =coord_setup.repl_force      ;

    size_t M_ves = VES_STRIDE;                 // Points per vesicle
    size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
    assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

    size_t N_trg=trg_coord.Dim()/COORD_DIM;
    proj_patch_param.ReInit(N_trg*        2);
    proj_coord      .ReInit(N_trg*COORD_DIM);
    repl_force      .ReInit(N_trg*COORD_DIM);
    is_extr_pt      .ReInit(N_trg          );
    Real_t r2repul_inv=(repul_dist_>0?std::pow(1.0/repul_dist_,2.0):0);
    Real_t& r_near=coord_setup.r_near;
    pvfmm::Vector<Real_t> min_dist_loc_(omp_p);
    min_dist_loc_.SetZero();
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Setup
      Real_t min_dist_loc=1e10;
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
        // Compute projection for near points
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          QuadraticPatch patch;
          { // create patch
            Real_t mesh[3*3*COORD_DIM];
            size_t k_=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
            QuadraticPatch::patch_mesh(mesh, sh_order_, k_, &S[0][i*M_ves*COORD_DIM]);
            patch=QuadraticPatch(&mesh[0],COORD_DIM);
          }
          { // Find nearest point on patch (first interpolation point)
            Real_t& x=proj_patch_param[trg_idx*2+0];
            Real_t& y=proj_patch_param[trg_idx*2+1];
            is_extr_pt[trg_idx]=
            patch.project(&  trg_coord[trg_idx*COORD_DIM],x,y);
            patch.eval(x,y,&proj_coord[trg_idx*COORD_DIM]);
            Real_t normal[COORD_DIM];
            { // Compute normal
              Real_t sgrad[2*COORD_DIM];
              patch.grad(x,y,sgrad);
              normal[0]=sgrad[1]*sgrad[5]-sgrad[2]*sgrad[4];
              normal[1]=sgrad[2]*sgrad[3]-sgrad[0]*sgrad[5];
              normal[2]=sgrad[0]*sgrad[4]-sgrad[1]*sgrad[3];
              Real_t r=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
              normal[0]/=r;
              normal[1]/=r;
              normal[2]/=r;
            }

            { // repulsion
              Real_t r=0;
              r+=(trg_coord[trg_idx*COORD_DIM+0]-proj_coord[trg_idx*COORD_DIM+0])*normal[0];
              r+=(trg_coord[trg_idx*COORD_DIM+1]-proj_coord[trg_idx*COORD_DIM+1])*normal[1];
              r+=(trg_coord[trg_idx*COORD_DIM+2]-proj_coord[trg_idx*COORD_DIM+2])*normal[2];
              min_dist_loc=std::min(min_dist_loc,r);

              Real_t f=0;
              if(r2repul_inv>0 && r<r_near) f=exp(-r*r*r2repul_inv)/r*(1.0-r/r_near);
              repl_force[trg_idx*COORD_DIM+0]=normal[0]*f;
              repl_force[trg_idx*COORD_DIM+1]=normal[1]*f;
              repl_force[trg_idx*COORD_DIM+2]=normal[2]*f;
            }
          }
        }
      }
      min_dist_loc_[tid]=min_dist_loc;
    }
    pvfmm::Profile::Toc();

    { // Print minimum distance between surfaces
      Real_t min_dist_loc=-1e10, min_dist;
      for(long i=0;i<omp_p;i++) min_dist_loc=std::max(min_dist_loc, -min_dist_loc_[i]);
      MPI_Allreduce(&min_dist_loc, &min_dist, 1, pvfmm::par::Mpi_datatype<Real_t>::value(), pvfmm::par::Mpi_datatype<Real_t>::max(), comm);
      INFO("Min-distance = "<<-min_dist<<'\n');
    }
  }

  { // repulsion
    PVFMMVec_t& repl_force=coord_setup.repl_force;
    VelocityScatter(repl_force);
    assert(repl_force.Dim()==T.Dim());
  }
  pvfmm::Profile::Enable(prof_state);
  pvfmm::Profile::Toc();
}

template<typename Real_t>
const typename NearSingular<Real_t>::PVFMMVec_t&  NearSingular<Real_t>::ForceRepul(){
  SetupCoordData();
  PVFMMVec_t& repl_force=coord_setup.repl_force;
  return repl_force;
}

template<typename Real_t>
bool NearSingular<Real_t>::CheckCollision(){
  SetupCoordData();
  bool collision=false;
  pvfmm::Vector<char>& is_extr_pt=coord_setup.is_extr_pt;
  for(long i=0;i<is_extr_pt.Dim();i++){
    if(!is_extr_pt[i]){
      collision=true;
      break;
    }
  }
  return collision;
}

template<typename Real_t>
void NearSingular<Real_t>::SubtractDirect(PVFMMVec_t& vel_fmm){
  if(update_direct){ // Compute vel_direct
    size_t omp_p=omp_get_max_threads();
    SetupCoordData();

    pvfmm::Profile::Tic("SubtractDirect",&comm,true);
    PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
    pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
    pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;
    vel_direct.ReInit(trg_coord.Dim()); vel_direct.SetZero();

    size_t M_ves = VES_STRIDE;                 // Points per vesicle
    size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
    assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Compute vel_direct for near points.
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
        PVFMMVec_t s_coord(     M_ves*COORD_DIM, const_cast<Real_t*>(&S[0]      [   i*M_ves*COORD_DIM]), false);
        PVFMMVec_t t_coord(trg_cnt[i]*COORD_DIM, &trg_coord [trg_dsp[i]*COORD_DIM], false);
        PVFMMVec_t t_veloc(trg_cnt[i]*COORD_DIM, &vel_direct[trg_dsp[i]*COORD_DIM], false);

        if(qforce_single){ // Subtract wrong near potential
          PVFMMVec_t qforce(M_ves*(COORD_DIM*1), const_cast<Real_t*>(&qforce_single[0][0]+M_ves*(COORD_DIM*1)*i), false);
          StokesKernel<Real_t>::Kernel().k_s2t->      ker_poten(&s_coord[0], M_ves, &qforce[0], 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
        }
        if(qforce_double){ // Subtract wrong near potential
          PVFMMVec_t qforce(M_ves*(COORD_DIM*2), const_cast<Real_t*>(&qforce_double[0][0]+M_ves*(COORD_DIM*2)*i), false);
          StokesKernel<Real_t>::Kernel().k_s2t->dbl_layer_poten(&s_coord[0], M_ves, &qforce[0], 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
        }
      }
    }
    VelocityScatter(vel_direct);

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++) if(trg_is_surf){ // Compute vel_direct for surface points.
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){ // loop over all vesicles
        size_t Ta=(vel_direct.Dim()/COORD_DIM)*(i+0)/N_ves;
        size_t Tb=(vel_direct.Dim()/COORD_DIM)*(i+1)/N_ves;

        PVFMMVec_t s_coord(  M_ves*COORD_DIM, const_cast<Real_t*>(&S[0]      [i*M_ves*COORD_DIM]), false);
        PVFMMVec_t t_coord((Tb-Ta)*COORD_DIM, &T         [     Ta*COORD_DIM], false);
        PVFMMVec_t t_veloc((Tb-Ta)*COORD_DIM, &vel_direct[     Ta*COORD_DIM], false);

        if(qforce_single){ // Subtract wrong near potential
          PVFMMVec_t qforce(M_ves*(COORD_DIM*1), const_cast<Real_t*>(&qforce_single[0][0]+M_ves*(COORD_DIM*1)*i), false);
          StokesKernel<Real_t>::Kernel().k_s2t->      ker_poten(&s_coord[0], M_ves, &qforce[0], 1, &t_coord[0], Tb-Ta, &t_veloc[0], NULL);
        }
        if(qforce_double){ // Subtract wrong near potential
          PVFMMVec_t qforce(M_ves*(COORD_DIM*2), const_cast<Real_t*>(&qforce_double[0][0]+M_ves*(COORD_DIM*2)*i), false);
          StokesKernel<Real_t>::Kernel().k_s2t->dbl_layer_poten(&s_coord[0], M_ves, &qforce[0], 1, &t_coord[0], Tb-Ta, &t_veloc[0], NULL);
        }
      }
    }

    update_direct=NearSingular::UpdateNone;
    pvfmm::Profile::Toc();
  }

  assert(vel_direct.Dim()==vel_fmm.Dim());
  #pragma omp parallel for
  for(size_t i=0;i<vel_direct.Dim();i++){
    vel_fmm[i]-=vel_direct[i];
  }
}

template<typename Real_t>
void NearSingular<Real_t>::VelocityScatter(PVFMMVec_t& trg_vel){
  PVFMMVec_t trg_vel_in;
  { // Initialize trg_vel_in <-- trg_vel; trg_vel <-- Zeros(trg_count*COORD_DIM);
    trg_vel_in.Swap(trg_vel);
    trg_vel.ReInit(T.Dim());
    trg_vel.SetZero();
  }
  { // Scatter trg_vel+=Scatter(trg_vel_in)
    pvfmm::Profile::Tic("ScatterTrg",&comm,true);
    pvfmm::Vector<size_t>& trg_scatter=coord_setup.near_trg_scatter;
    pvfmm::Vector<size_t>&   trg_pt_id=coord_setup.near_trg_pt_id;
    pvfmm::par::ScatterForward(trg_vel_in, trg_scatter, comm);

    size_t trg_id_offset;
    { // Get trg_id_offset
      long long disp=0;
      long long size=trg_vel.Dim()/COORD_DIM;
      MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
      trg_id_offset=disp-size;
    }

    size_t omp_p=omp_get_max_threads();
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){
      size_t a=((tid+0)*trg_pt_id.Dim())/omp_p;
      size_t b=((tid+1)*trg_pt_id.Dim())/omp_p;
      while(a>0 && a<trg_pt_id.Dim() && trg_pt_id[a-1]==trg_pt_id[a]) a++;
      while(b>0 && b<trg_pt_id.Dim() && trg_pt_id[b-1]==trg_pt_id[b]) b++;
      for(size_t i=a;i<b;i++){
        size_t pt_id=trg_pt_id[i]-trg_id_offset;
        trg_vel[pt_id*COORD_DIM+0]+=trg_vel_in[i*COORD_DIM+0];;
        trg_vel[pt_id*COORD_DIM+1]+=trg_vel_in[i*COORD_DIM+1];;
        trg_vel[pt_id*COORD_DIM+2]+=trg_vel_in[i*COORD_DIM+2];;
      }
    }
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
    node_data.max_pts=coord_setup.near_trg_coord.Dim();

    { // Set node_data.pt_coord, node_data.pt_value
      PVFMMVec_t&           near_trg_coord=coord_setup.near_trg_coord;
      pvfmm::Vector<size_t>&  near_trg_cnt=coord_setup.  near_trg_cnt;

      Real_t scale_x, shift_x[COORD_DIM];
      {
        shift_x[0]=coord_setup.bbox[0];
        shift_x[1]=coord_setup.bbox[1];
        shift_x[2]=coord_setup.bbox[2];
        scale_x=coord_setup.bbox[3];
      }

      node_data.pt_coord=near_trg_coord;
      for(size_t i=0;i<node_data.pt_coord.Dim()/COORD_DIM;i++){
        for(size_t k=0;k<COORD_DIM;k++){
          Real_t c=node_data.pt_coord[i*COORD_DIM+k];
          c=c*scale_x+shift_x[k];
          while(c< 0.0) c+=1.0;
          while(c>=1.0) c-=1.0;
          node_data.pt_coord[i*COORD_DIM+k]=c;
        }
      }

      std::vector<Real_t> value;
      size_t ves_id_offset;
      { // Get ves_id_offset
        size_t M_ves = VES_STRIDE;                 // Points per vesicle
        size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
        assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

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
    node_data.max_pts=T.Dim();

    { // Set node_data.pt_coord, node_data.pt_value
      Real_t scale_x, shift_x[COORD_DIM];
      {
        shift_x[0]=coord_setup.bbox[0];
        shift_x[1]=coord_setup.bbox[1];
        shift_x[2]=coord_setup.bbox[2];
        scale_x=coord_setup.bbox[3];
      }

      node_data.pt_coord=T;
      for(size_t i=0;i<node_data.pt_coord.Dim()/COORD_DIM;i++){
        for(size_t k=0;k<COORD_DIM;k++){
          Real_t c=node_data.pt_coord[i*COORD_DIM+k];
          c=c*scale_x+shift_x[k];
          while(c< 0.0) c+=1.0;
          while(c>=1.0) c-=1.0;
          node_data.pt_coord[i*COORD_DIM+k]=c;
        }
      }
      node_data.pt_value=trg_vel;
    }

    pt_tree.Initialize(&node_data);
    pt_tree.Write2File("trg_vel");
    pvfmm::Profile::Toc();
  }
}



template <class Real_t>
static inline Real_t InterPoints(int i, int deg){
  return (i?1.0+(i-1.0)/(deg-2):0.0);  // Clustered points
  //return i;                            // Equi-spaced points
}

template <class Real_t>
static inline Real_t InterPoly(Real_t x, int j, int deg){
  Real_t y=1.0;
  Real_t x0=InterPoints<Real_t>(j,deg);
  for(size_t k=0;k<deg;k++){ // Lagrange polynomial
    Real_t xk=InterPoints<Real_t>(k,deg);
    if(j!=k) y*=(x-xk)/(x0-xk);
  }
  //y=pow(x,j)                       // Polynomial basis
  //y=pow(1.0/(1.0+x),j+1)           // Laurent polynomial
  return y;
}

template<typename Real_t>
const typename NearSingular<Real_t>::PVFMMVec_t& NearSingular<Real_t>::operator()(bool update){
  if(!update || !update_interp){
    return vel_interp;
  }
  update_interp=NearSingular::UpdateNone;
  static Real_t eps_sqrt=-1;
  if(eps_sqrt<0){
    #pragma omp critical
    if(eps_sqrt<0){
      eps_sqrt=1.0;
      while(eps_sqrt+(Real_t)1.0>1.0) eps_sqrt*=0.5;
      eps_sqrt=sqrt(eps_sqrt);
    }
  }

  pvfmm::Profile::Tic("NearInteraction",&comm,true);
  bool prof_state=pvfmm::Profile::Enable(false);
  size_t omp_p=omp_get_max_threads();
  SetupCoordData();
  assert(S_vel);

  Real_t&                  r_near=coord_setup.        r_near;
  PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
  pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
  pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;

  pvfmm::Vector<char>& is_extr_pt=coord_setup.is_extr_pt;
  PVFMMVec_t& proj_patch_param=coord_setup.proj_patch_param;
  PVFMMVec_t& proj_coord      =coord_setup.proj_coord      ;

  size_t M_ves = VES_STRIDE;                 // Points per vesicle
  size_t N_ves = S->Dim()/(M_ves*COORD_DIM); // Number of vesicles
  assert(N_ves*(M_ves*COORD_DIM) == S->Dim());

  vel_interp.ReInit(trg_coord.Dim());
  pvfmm::Profile::Tic("VelocInterp",&comm,true);
  #pragma omp parallel for
  for(size_t tid=0;tid<omp_p;tid++){ // Compute vel_interp.
    PVFMMVec_t interp_coord;
    PVFMMVec_t interp_veloc;
    PVFMMVec_t patch_veloc;
    PVFMMVec_t interp_x;

    size_t a=((tid+0)*N_ves)/omp_p;
    size_t b=((tid+1)*N_ves)/omp_p;
    for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
      PVFMMVec_t s_coord(M_ves*COORD_DIM, const_cast<Real_t*>(&S[0][i*M_ves*COORD_DIM]), false);
      { // Resize interp_coord, interp_veloc, patch_veloc, interp_x; interp_veloc[:]=0
        interp_coord.Resize(trg_cnt[i]*(INTERP_DEG-1)*COORD_DIM);
        interp_veloc.Resize(trg_cnt[i]*(INTERP_DEG-1)*COORD_DIM);
        patch_veloc .Resize(trg_cnt[i]               *COORD_DIM);
        interp_x    .Resize(trg_cnt[i]                         );
        interp_veloc.SetZero();
      }
      { // Set interp_x, interp_coord
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          Real_t interp_coord0[COORD_DIM]={proj_coord[trg_idx*COORD_DIM+0],
                                           proj_coord[trg_idx*COORD_DIM+1],
                                           proj_coord[trg_idx*COORD_DIM+2]};
          Real_t dR[COORD_DIM]={trg_coord[trg_idx*COORD_DIM+0]-interp_coord0[0],
                                trg_coord[trg_idx*COORD_DIM+1]-interp_coord0[1],
                                trg_coord[trg_idx*COORD_DIM+2]-interp_coord0[2]};
          Real_t dR_norm=sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
          if(!is_extr_pt[trg_idx] && trg_is_surf) dR_norm=-dR_norm; // always use exterior normal
          Real_t OOdR=1.0/dR_norm;
          if(fabs(dR_norm)<eps_sqrt){
            dR_norm=0;
            OOdR=0;
          }

          interp_x[j]=dR_norm/r_near;
          for(size_t l=0;l<INTERP_DEG-1;l++){
            Real_t x=InterPoints<Real_t>(l+1, INTERP_DEG);
            interp_coord[(l+j*(INTERP_DEG-1))*COORD_DIM+0]=interp_coord0[0]+dR[0]*OOdR*r_near*x;
            interp_coord[(l+j*(INTERP_DEG-1))*COORD_DIM+1]=interp_coord0[1]+dR[1]*OOdR*r_near*x;
            interp_coord[(l+j*(INTERP_DEG-1))*COORD_DIM+2]=interp_coord0[2]+dR[2]*OOdR*r_near*x;
          }
        }
      }
      if(interp_veloc.Dim()){ // Set interp_veloc
        if(qforce_single){
          StokesKernel<Real_t>::Kernel().k_s2t->ker_poten(&s_coord[0], M_ves, const_cast<Real_t*>(&qforce_single[0][0]+M_ves*(COORD_DIM*1)*i), 1, &interp_coord[0], interp_coord.Dim()/COORD_DIM, &interp_veloc[0], NULL);
        }
        if(qforce_double){
          StokesKernel<Real_t>::Kernel().k_s2t->dbl_layer_poten(&s_coord[0], M_ves, const_cast<Real_t*>(&qforce_double[0][0]+M_ves*(COORD_DIM*2)*i), 1, &interp_coord[0], interp_coord.Dim()/COORD_DIM, &interp_veloc[0], NULL);
        }
      }
      { // Set patch_veloc
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          QuadraticPatch patch;
          { // create patch
            Real_t mesh[3*3*COORD_DIM];
            { // compute mesh
              size_t k_=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
              QuadraticPatch::patch_mesh(mesh   , sh_order_, k_, &S_vel       [0][i*M_ves*COORD_DIM]);
            }
            if(force_double){ // add contribution from force_double to mesh
              Real_t mesh_fd[3*3*COORD_DIM];
              size_t k_=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
              QuadraticPatch::patch_mesh(mesh_fd, sh_order_, k_, &force_double[0][i*M_ves*COORD_DIM]);
              Real_t scal=0.5; if(!is_extr_pt[trg_idx] && !trg_is_surf) scal=-0.5;
              for(size_t k=0;k<3*3*COORD_DIM;k++) mesh[k]+=scal*mesh_fd[k];
            }
            patch=QuadraticPatch(&mesh[0],COORD_DIM);
          }
          { // compute first interpolation point
            Real_t& x=proj_patch_param[trg_idx*2+0];
            Real_t& y=proj_patch_param[trg_idx*2+1];
            patch.eval(x,y,&patch_veloc[j*COORD_DIM]);
          }
        }
      }

      { // Interpolate
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          Real_t* veloc_interp_=&vel_interp[trg_idx*COORD_DIM];
          if(interp_x[j]==0){
            for(size_t k=0;k<COORD_DIM;k++){
              veloc_interp_[k]=patch_veloc[j*COORD_DIM+k];
            }
          }else{ // Interpolate and find target velocity veloc_interp_
            pvfmm::Matrix<Real_t> y(INTERP_DEG,1);
            pvfmm::Matrix<Real_t> x(INTERP_DEG,1);
            Real_t x_=interp_x[j];
            for(size_t l=0;l<INTERP_DEG;l++){
              x[l][0]=InterPoly(x_,l,INTERP_DEG);
            }
            for(size_t k=0;k<COORD_DIM;k++){
              y[0][0]=patch_veloc[j*COORD_DIM+k];
              for(size_t l=0;l<INTERP_DEG-1;l++){
                y[l+1][0]=interp_veloc[(l+j*(INTERP_DEG-1))*COORD_DIM+k];
              }
              #if 0
              static pvfmm::Matrix<Real_t> M;
              if(M.Dim(0)*M.Dim(1)==0){ // matrix for computing interpolation coefficients
                #pragma omp critical
                if(M.Dim(0)*M.Dim(1)==0){
                  pvfmm::Matrix<Real_t> M_(INTERP_DEG,INTERP_DEG);
                  assert(InterPoints<Real_t>(0, INTERP_DEG)==0);
                  for(size_t i=0;i<INTERP_DEG;i++){
                    Real_t x=InterPoints<Real_t>(i, INTERP_DEG);
                    for(size_t j=0;j<INTERP_DEG;j++){
                      M_[i][j]=InterPoly(x,j,INTERP_DEG);
                    }
                  }
                  M=M_.pinv();
                }
              }
              pvfmm::Matrix<Real_t> coeff(INTERP_DEG,1);
              coeff=M*y;
              #else
              pvfmm::Matrix<Real_t>& coeff=y;
              #endif
              veloc_interp_[k]=0;
              for(size_t l=0;l<INTERP_DEG;l++) veloc_interp_[k]+=y[l][0]*x[l][0];
            }
          }
        }
      }
    }
  }
  pvfmm::Profile::Toc();

  pvfmm::Profile::Tic("Scatter",&comm,true);
  VelocityScatter(vel_interp);
  assert(vel_interp.Dim()==T.Dim());
  pvfmm::Profile::Toc();

  pvfmm::Profile::Enable(prof_state);
  pvfmm::Profile::Toc();

  return vel_interp;
}



template<typename Real_t>
NearSingular<Real_t>::QuadraticPatch::QuadraticPatch(Real_t* x, int dof_){
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

template<typename Real_t>
void NearSingular<Real_t>::QuadraticPatch::eval(Real_t x, Real_t y, Real_t* val){
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

template<typename Real_t>
int NearSingular<Real_t>::QuadraticPatch::project(Real_t* t_coord_j, Real_t& x, Real_t&y){ // Find nearest point on patch
  static Real_t eps=-1;
  if(eps<0){
    #pragma omp critical
    if(eps<0){
      eps=1.0;
      while(eps+(Real_t)1.0>1.0) eps*=0.5;
    }
  }

  assert(dof==COORD_DIM);
  Real_t sc[COORD_DIM];
  x=0; y=0;
  eval(x,y,sc);
  while(1){
    Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                          t_coord_j[1]-sc[1],
                          t_coord_j[2]-sc[2]};
    Real_t dR2=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];

    Real_t sgrad[2*COORD_DIM];
    grad(x,y,sgrad);

    { // Break when |dR-(dR.n)n| < 1e-8
      Real_t n[COORD_DIM];
      n[0]=sgrad[1]*sgrad[3+2]-sgrad[2]*sgrad[3+1];
      n[1]=sgrad[2]*sgrad[3+0]-sgrad[0]*sgrad[3+2];
      n[2]=sgrad[0]*sgrad[3+1]-sgrad[1]*sgrad[3+0];
      Real_t n_norm=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      n[0]/=n_norm;
      n[1]/=n_norm;
      n[2]/=n_norm;

      Real_t dRn=dR[0]*n[0]+dR[1]*n[1]+dR[2]*n[2];
      dR[0]-=dRn*n[0];
      dR[1]-=dRn*n[1];
      dR[2]-=dRn*n[2];

      Real_t norm=sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
      if(norm<1e-8) break;
    }
    Real_t dxdx=sgrad[0]*sgrad[0]+sgrad[1]*sgrad[1]+sgrad[2]*sgrad[2];
    Real_t dydy=sgrad[3]*sgrad[3]+sgrad[4]*sgrad[4]+sgrad[5]*sgrad[5];
    Real_t dxdy=sgrad[0]*sgrad[3]+sgrad[1]*sgrad[4]+sgrad[2]*sgrad[5];
    Real_t dxdR=sgrad[0]*   dR[0]+sgrad[1]*   dR[1]+sgrad[2]*   dR[2];
    Real_t dydR=sgrad[3]*   dR[0]+sgrad[4]*   dR[1]+sgrad[5]*   dR[2];

    Real_t dx, dy;
    { // Set dx, dy
      dx=(dydy*dxdR-dxdy*dydR)/(dxdx*dydy-dxdy*dxdy);
      dy=(dxdx*dydR-dxdy*dxdR)/(dxdx*dydy-dxdy*dxdy);
    }
    if(1){ // Check for cases which should not happen and break;
      if((x<=-1.2 && dx<0.0) ||
         (x>= 1.2 && dx>0.0) ||
         (y<=-1.2 && dy<0.0) ||
         (y>= 1.2 && dy>0.0)){
        break;
      }
    }
    if(1){ // if x+dx or y+dy are outsize [-1.2,1.2]
      if(x>-1.2 && x+dx<-1.2){
        Real_t s=(-1.2-x)/dx;
        dx*=s; dy*=s;
      }
      if(x< 1.2 && x+dx> 1.2){
        Real_t s=( 1.2-x)/dx;
        dx*=s; dy*=s;
      }
      if(y>-1.2 && y+dy<-1.2){
        Real_t s=(-1.2-y)/dy;
        dx*=s; dy*=s;
      }
      if(y< 1.2 && y+dy> 1.2){
        Real_t s=( 1.2-y)/dy;
        dx*=s; dy*=s;
      }
    }

    while(1){ // increment x,y
      { // Account for curvature.
        Real_t dR2_[3]={dR2,0,0};
        { // x+dx*0.5,y+dy*0.5
          Real_t sc[COORD_DIM];
          eval(x+dx*0.5,y+dy*0.5,sc);
          Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                                t_coord_j[1]-sc[1],
                                t_coord_j[2]-sc[2]};
          dR2_[1]=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
        }
        { // x+dx*1.0,y+dy*1.0
          Real_t sc[COORD_DIM];
          eval(x+dx*1.0,y+dy*1.0,sc);
          Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                                t_coord_j[1]-sc[1],
                                t_coord_j[2]-sc[2]};
          dR2_[2]=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
        }

        Real_t coeff[3];
        coeff[0]=                         dR2_[0]    ;
        coeff[1]=-dR2_[2]    +dR2_[1]*4.0-dR2_[0]*3.0;
        coeff[2]= dR2_[2]*2.0-dR2_[1]*4.0+dR2_[0]*2.0;
        if(coeff[2]>0){
          Real_t s=-0.5*coeff[1]/coeff[2];
          if(s>0.0 && s<1.0){
            dx*=s; dy*=s;
          }
        }
      }
      { // check if dx, dy reduce dR2
        Real_t dR2_;
        eval(x+dx,y+dy,sc);
        Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                              t_coord_j[1]-sc[1],
                              t_coord_j[2]-sc[2]};
        dR2_=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
        if(dR2_!=dR2_) assert(false); // Check NaN
        if(dR2_<dR2){ x+=dx; y+=dy; break;}
        if(dx*dx*dxdx+dy*dy*dydy<64*eps) break;
        else {dx*=0.5; dy*=0.5;}
      }
    }
    if(dx*dx*dxdx+dy*dy*dydy<64*eps) break;
  }

  { // Determine direction of point (exterior or interior)
    Real_t dR[COORD_DIM]={t_coord_j[0]-sc[0],
                          t_coord_j[1]-sc[1],
                          t_coord_j[2]-sc[2]};
    Real_t sgrad[2*COORD_DIM];
    grad(x,y,sgrad);

    Real_t direc=0;
    direc+=dR[0]*sgrad[1]*sgrad[3+2];
    direc+=dR[1]*sgrad[2]*sgrad[3+0];
    direc+=dR[2]*sgrad[0]*sgrad[3+1];
    direc-=dR[0]*sgrad[2]*sgrad[3+1];
    direc-=dR[1]*sgrad[0]*sgrad[3+2];
    direc-=dR[2]*sgrad[1]*sgrad[3+0];
    return (direc>0);
  }
}

template<typename Real_t>
void NearSingular<Real_t>::QuadraticPatch::patch_mesh(Real_t* patch_value_, size_t sh_order, size_t k_, const Real_t* s_value){
  static pvfmm::Matrix<Real_t> M_patch_interp[SHMAXDEG];
  { // Compute interpolation matrix for the pole
    static std::vector<bool> flag(SHMAXDEG,false);
    assert(sh_order<SHMAXDEG);
    if(!flag[sh_order]){
      #pragma omp critical (NEARSINGULAR_PATCH_INTERP)
      if(!flag[sh_order]){
        size_t k1_max=2*sh_order;
        pvfmm::Matrix<Real_t>& M=M_patch_interp[sh_order];
        M.ReInit(8,k1_max);
        M.SetZero();
        for(size_t k0=0;k0<8;k0++)
        for(size_t k1=0;k1<k1_max;k1++)
        for(size_t k2=0;k2<=sh_order;k2++){
          M[k0][k1]+=(!k2 || k2==sh_order?1.0:2.0)*sin(2.0*M_PI*k2*(k0*1.0/8.0))*sin(2.0*M_PI*k2*(k1*1.0/k1_max))*1.0/k1_max;
          M[k0][k1]+=(!k2 || k2==sh_order?1.0:2.0)*cos(2.0*M_PI*k2*(k0*1.0/8.0))*cos(2.0*M_PI*k2*(k1*1.0/k1_max))*1.0/k1_max;
        }
        flag[sh_order]=true;
      }
    }
  }

  size_t k0_max=1+sh_order;
  size_t k1_max=2*sh_order;
  if(k_==0){ // Create patch
    Real_t tmp[COORD_DIM*8];
    pvfmm::Matrix<Real_t> value(8,COORD_DIM,tmp,false);
    { // Set value
      int k=2+0*k1_max;
      pvfmm::Matrix<Real_t> X(k1_max,COORD_DIM,(Real_t*)&s_value[k*COORD_DIM],false);
      pvfmm::Matrix<Real_t>::GEMM(value,M_patch_interp[sh_order],X);
    }

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+0)*COORD_DIM+i]=value[0][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+1)*COORD_DIM+i]=value[7][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+2)*COORD_DIM+i]=value[6][i];

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+0)*COORD_DIM+i]=value[1][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+1)*COORD_DIM+i]=s_value[0*COORD_DIM+i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+2)*COORD_DIM+i]=value[5][i];

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+0)*COORD_DIM+i]=value[2][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+1)*COORD_DIM+i]=value[3][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+2)*COORD_DIM+i]=value[4][i];

  }else if(k_==1){
    Real_t tmp[COORD_DIM*8];
    pvfmm::Matrix<Real_t> value(8,COORD_DIM,tmp,false);
    { // Set value
      int k=2+(k0_max-1)*k1_max;
      pvfmm::Matrix<Real_t> X(k1_max,COORD_DIM,(Real_t*)&s_value[k*COORD_DIM],false);
      pvfmm::Matrix<Real_t>::GEMM(value,M_patch_interp[sh_order],X);
    }

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+0)*COORD_DIM+i]=value[0][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+1)*COORD_DIM+i]=value[1][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+2)*COORD_DIM+i]=value[2][i];

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+0)*COORD_DIM+i]=value[7][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+1)*COORD_DIM+i]=s_value[1*COORD_DIM+i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+2)*COORD_DIM+i]=value[3][i];

    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+0)*COORD_DIM+i]=value[6][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+1)*COORD_DIM+i]=value[5][i];
    for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+2)*COORD_DIM+i]=value[4][i];

  }else{
    int k0_, k1_; // mesh coordinates for nearest point
    { // Find nearest point on mesh
      k0_=(k_-2)/k1_max;
      k1_=(k_-2)%k1_max;
    }

    { // (0,0)
      int k=0; if(k0_>0) k=2+((k1_+k1_max-1)%k1_max)+(k0_-1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+0)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (0,1)
      int k=0; if(k0_>0) k=2+  k1_                  +(k0_-1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+1)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (0,2)
      int k=0; if(k0_>0) k=2+((k1_       +1)%k1_max)+(k0_-1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(0*3+2)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }

    { // (1,0)
      int k=2+((k1_+k1_max-1)%k1_max)+k0_*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+0)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (1,1)
      int k=2+  k1_                  +k0_*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+1)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (1,2)
      int k=2+((k1_       +1)%k1_max)+k0_*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(1*3+2)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }

    { // (2,0)
      int k=1; if(k0_<k0_max-1) k=2+((k1_+k1_max-1)%k1_max)+(k0_+1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+0)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (2,1)
      int k=1; if(k0_<k0_max-1) k=2+  k1_                  +(k0_+1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+1)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }
    { // (2,2)
      int k=1; if(k0_<k0_max-1) k=2+((k1_       +1)%k1_max)+(k0_+1)*k1_max;
      for(size_t i=0;i<COORD_DIM;i++) patch_value_[(2*3+2)*COORD_DIM+i]=s_value[k*COORD_DIM+i];
    }

  }

  if (k_==0 || k_==1) { // Fix patch at poles
    static pvfmm::Matrix<Real_t> M;
    static bool compute_M = true;
    if(compute_M){
      #pragma omp critical
      if(compute_M){
        pvfmm::Matrix<Real_t> M0(3*3,3*3), M1(3*3,3*3);
        for(int i=0;i<3;i++){ // Set M0, M1
          for(int j=0;j<3;j++){
            Real_t x0 = i - 1, y0 = j - 1;
            Real_t x1 = x0, y1 = y0;
            if(x1 && y1){
              x1=x1/sqrt(2);
              y1=y1/sqrt(2);
            }
            for(int m=0;m<3;m++){
              for(int n=0;n<3;n++){
                M0[i*3+j][m*3+n]=pvfmm::pow<Real_t>(x0,m)*pvfmm::pow<Real_t>(y0,n);
                M1[i*3+j][m*3+n]=pvfmm::pow<Real_t>(x1,m)*pvfmm::pow<Real_t>(y1,n);
              }
            }
          }
        }
        M = M0 * M1.pinv();
        compute_M = false;
      }
      #pragma omp flush(compute_M)
    }
    pvfmm::Matrix<Real_t> v(3*3, COORD_DIM, patch_value_, false);
    v = M * v;
  }
}

template<typename Real_t>
void NearSingular<Real_t>::QuadraticPatch::grad(Real_t x, Real_t y, Real_t* val){
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

#undef VES_STRIDE
