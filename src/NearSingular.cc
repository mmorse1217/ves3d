#ifdef HAVE_PVFMM

#include <omp.h>
#include <iostream>

#include <ompUtils.h>
#include <parUtils.h>
#include <profile.hpp>
#include <mortonid.hpp>
#include <legendre_rule.hpp>
#include <mpi_tree.hpp> // Only for vis

template<typename Surf_t>
NearSingular<Surf_t>::NearSingular(Real_t box_size, MPI_Comm c){
  box_size_=box_size;
  comm=c;
  S=NULL;
  qforce_single=NULL;
  qforce_double=NULL;
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


template<typename Surf_t>
void NearSingular<Surf_t>::SetSrcCoord(const Surf_t& S_){
  update_direct=update_direct | NearSingular::UpdateSrcCoord;
  update_interp=update_interp | NearSingular::UpdateSrcCoord;
  update_setup =update_setup  | NearSingular::UpdateSrcCoord;
  S=&S_;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetSurfaceVel(const Vec_t& S_vel_){
  update_interp=update_interp | NearSingular::UpdateSurfaceVel;
  S_vel=&S_vel_;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetDensitySL(const PVFMMVec_t* qforce_single_){
  update_direct=update_direct | NearSingular::UpdateDensitySL;
  update_interp=update_interp | NearSingular::UpdateDensitySL;
  qforce_single=qforce_single_;
}

template<typename Surf_t>
void NearSingular<Surf_t>::SetDensityDL(const PVFMMVec_t* qforce_double_, const Vec_t* force_double_){
  update_direct=update_direct | NearSingular::UpdateDensityDL;
  update_interp=update_interp | NearSingular::UpdateDensityDL;
  qforce_double=qforce_double_;
  force_double=force_double_;
}


template<typename Surf_t>
void NearSingular<Surf_t>::SetTrgCoord(Real_t* trg_coord, size_t N){ // TODO: change to const Real_t*
  update_direct=update_direct | NearSingular::UpdateTrgCoord;
  update_interp=update_interp | NearSingular::UpdateTrgCoord;
  update_setup =update_setup  | NearSingular::UpdateTrgCoord;
  T.ReInit(N*COORD_DIM,trg_coord);
}



template<typename Surf_t>
void NearSingular<Surf_t>::SetupCoordData(){
  assert(S);
  Real_t near=2.0/sqrt((Real_t)S->getShOrder()); // TODO: some function of sh_order and accuracy
  if(!(update_setup & (NearSingular::UpdateSrcCoord | NearSingular::UpdateTrgCoord))) return;
  update_setup=update_setup & ~(NearSingular::UpdateSrcCoord | NearSingular::UpdateTrgCoord);

  int np, rank;
  MPI_Comm_size(comm,&np);
  MPI_Comm_rank(comm,&rank);
  size_t omp_p=omp_get_max_threads();

  assert(S);
  struct{
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
    Real_t& r_near=coord_setup.r_near;
    { // Set tree pt data
      pvfmm::Profile::Tic("PtData",&comm,true);
      Real_t* bbox=coord_setup.bbox;
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
      if(box_size_>0 && 2*r_near+r_ves>box_size_){
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
          PVFMMBoundingBox(      N_ves*M_ves, &pt_coord[0], &s0, x0, comm);
          PVFMMBoundingBox(T.Dim()/COORD_DIM, &       T[0], &s1, x1, comm);

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
              c[k]=pt_coord[i*COORD_DIM+k]*scale_x+shift_x[k];
              while(c[k]< 0.0) c[k]+=1.0;
              while(c[k]>=1.0) c[k]-=1.0;
            }
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

  { // find vesicle, near target pairs
    pvfmm::Profile::Tic("TrgNear",&comm,true);

    PVFMMVec_t&            near_trg_coord  =coord_setup.near_trg_coord  ;
    pvfmm::Vector<size_t>& near_trg_cnt    =coord_setup.near_trg_cnt    ;
    pvfmm::Vector<size_t>& near_trg_dsp    =coord_setup.near_trg_dsp    ;
    pvfmm::Vector<size_t>& near_trg_scatter=coord_setup.near_trg_scatter;
    pvfmm::Vector<size_t>& near_trg_pt_id  =coord_setup.near_trg_pt_id  ;
    pvfmm::Vector<size_t>& near_ves_pt_id  =coord_setup.near_ves_pt_id  ;
    // near_ves_pt_id: The neares vesicle point to the target point
    // TODO: This is not used right now, and the neares vesicle point is
    // compute again when constructing quadratic surface patch
    // // I think this is done now.

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

    { // Find near trg, ves points
      pvfmm::Profile::Tic("NearPair",&comm,true);
      std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_(omp_p); // (loc_trg_id, S_let.pt_id)
      #pragma omp parallel for
      for(size_t tid=0;tid<omp_p;tid++){
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_[tid];
        size_t tree_depth; Real_t r2_near;
        { // Set tree_depth, r_near
          tree_depth=S_let.mid[0].GetDepth();
          r2_near=coord_setup.r_near;
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
      const Vec_t& x=S->getPosition();
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

      const Vec_t& x=S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle

      for(size_t i=0;i<N_ves;i++){ // loop over all vesicles
        const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
        const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
        const Real_t* zk=x.getSubN_begin(i)+2*M_ves;

        for(size_t j=0;j<trg_cnt[i];j++){ // loop over near tagets
          size_t trg_idx=trg_dsp[i]+j;
          size_t src_idx=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
          Real_t ves_c[COORD_DIM]={xk[src_idx],yk[src_idx],zk[src_idx]};
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
}

template<typename Surf_t>
void NearSingular<Surf_t>::SubtractDirect(PVFMMVec_t& vel_fmm){
  if(update_direct){ // Compute vel_direct
    size_t omp_p=omp_get_max_threads();
    SetupCoordData();

    pvfmm::Profile::Tic("SubtractDirect",&comm,true);
    Real_t&                  r_near=coord_setup.        r_near;
    PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
    pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
    pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;
    vel_direct.ReInit(trg_coord.Dim()); vel_direct.SetZero();

    const Vec_t& x=S->getPosition();
    size_t N_ves = x.getNumSubs(); // Number of vesicles
    size_t M_ves = x.getStride(); // Points per vesicle
    int k0_max=x.getGridDim().first;
    int k1_max=x.getGridDim().second;
    assert(M_ves==k0_max*k1_max);

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Compute vel_direct.
      PVFMMVec_t s_coord(M_ves*COORD_DIM);

      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++){ // loop over all vesicles
        if(qforce_single || qforce_double){ // Set s_coord
          // read each component of x
          const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* zk=x.getSubN_begin(i)+2*M_ves;
          for(size_t j=0;j<M_ves;j++){
            s_coord[j*COORD_DIM+0]=xk[j];
            s_coord[j*COORD_DIM+1]=yk[j];
            s_coord[j*COORD_DIM+2]=zk[j];
          }
        }

        PVFMMVec_t t_coord(trg_cnt[i]*COORD_DIM, &trg_coord [trg_dsp[i]*COORD_DIM], false);
        PVFMMVec_t t_veloc(trg_cnt[i]*COORD_DIM, &vel_direct[trg_dsp[i]*COORD_DIM], false);
        t_veloc.SetZero();

        if(qforce_single){ // Subtract wrong near potential
          stokes_sl(&s_coord[0], M_ves, &qforce_single[0][0]+M_ves*(COORD_DIM*1)*i, 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
        }
        if(qforce_double){ // Subtract wrong near potential
          stokes_dl(&s_coord[0], M_ves, &qforce_double[0][0]+M_ves*(COORD_DIM*2)*i, 1, &t_coord[0], trg_cnt[i], &t_veloc[0], NULL);
        }
      }
    }

    VelocityScatter(vel_direct);
    update_direct=NearSingular::UpdateNone;
    pvfmm::Profile::Toc();
  }

  assert(vel_direct.Dim()==vel_fmm.Dim());
  #pragma omp parallel for
  for(size_t i=0;i<vel_direct.Dim();i++){
    vel_fmm[i]-=vel_direct[i];
  }
}

template<typename Surf_t>
void NearSingular<Surf_t>::VelocityScatter(PVFMMVec_t& trg_vel){
  if(coord_setup.near_trg_pt_id.Dim()){ // Scatter trg velocity
    pvfmm::Profile::Tic("ScatterTrg",&comm,true);

    size_t omp_p=omp_get_max_threads();
    pvfmm::Vector<size_t>& trg_scatter=coord_setup.near_trg_scatter;
    pvfmm::Vector<size_t>&   trg_pt_id=coord_setup.near_trg_pt_id;

    pvfmm::par::ScatterForward(trg_vel, trg_scatter, comm);

    size_t trg_count=T.Dim()/COORD_DIM;
    PVFMMVec_t trg_vel_final(trg_count*COORD_DIM);
    trg_vel_final.SetZero();

    size_t trg_id_offset;//=trg_pt_id[0];
    { // Get trg_id_offset
      long long disp=0;
      long long size=trg_count;
      MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
      trg_id_offset=disp-size;
    }

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
  }else{
    size_t trg_count=T.Dim()/COORD_DIM;
    trg_vel.ReInit(trg_count*COORD_DIM);
    trg_vel.SetZero();
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

      const Vec_t& x=S->getPosition();
      size_t N_ves = x.getNumSubs(); // Number of vesicles
      size_t M_ves = x.getStride(); // Points per vesicle
      assert(M_ves==x.getGridDim().first*x.getGridDim().second);

      std::vector<Real_t> value;
      size_t ves_id_offset;
      { // Scatter trg points to vesicles
        const Vec_t& x=S->getPosition();
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

  int project(Real_t* t_coord_j, Real_t& x, Real_t&y){ // Find nearest point on patch
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

      Real_t dxdx=sgrad[0]*sgrad[0]+sgrad[1]*sgrad[1]+sgrad[2]*sgrad[2];
      Real_t dydy=sgrad[3]*sgrad[3]+sgrad[4]*sgrad[4]+sgrad[5]*sgrad[5];
      Real_t dxdy=sgrad[0]*sgrad[3]+sgrad[1]*sgrad[4]+sgrad[2]*sgrad[5];
      Real_t dxdR=sgrad[0]*   dR[0]+sgrad[1]*   dR[1]+sgrad[2]*   dR[2];
      Real_t dydR=sgrad[3]*   dR[0]+sgrad[4]*   dR[1]+sgrad[5]*   dR[2];
      if(dxdR*dxdR/dxdx+dydR*dydR/dydy < dR2*(dxdx+dydy)*1e-6) break;

      Real_t dx, dy;
      if(dxdy){
        dx=(dxdR/dxdy-dydR/dydy)/(dxdx/dxdy-dxdy/dydy);
        dy=(dydR/dxdy-dxdR/dxdx)/(dydy/dxdy-dxdy/dxdx);
      }else{
        dx=dxdR/dxdx;
        dy=dydR/dydy;
      }
      { // Check for cases which should not happen and break;
        if((x<=-1.0 && dx<0.0) ||
           (x>= 1.0 && dx>0.0) ||
           (y<=-1.0 && dy<0.0) ||
           (y>= 1.0 && dy>0.0)){
          break;
        }
      }
      { // if x+dx or y+dy are outsize [-1,1]
        if(x>-1.0 && x+dx<-1.0){
          Real_t s=(-1.0-x)/dx;
          dx*=s; dy*=s;
        }
        if(x< 1.0 && x+dx> 1.0){
          Real_t s=( 1.0-x)/dx;
          dx*=s; dy*=s;
        }
        if(y>-1.0 && y+dy<-1.0){
          Real_t s=(-1.0-y)/dy;
          dx*=s; dy*=s;
        }
        if(y< 1.0 && y+dy> 1.0){
          Real_t s=( 1.0-y)/dy;
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
      return (direc<0);
    }
  }

  private:

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



template <class Real_t>
static void patch_mesh(Real_t* patch_value_, size_t sph_order, size_t k0_, size_t k1_,
    const Real_t* sx_coord, const Real_t* sy_coord, const Real_t* sz_coord, const Real_t* pole_coord, const Real_t* tcoord,
    const Real_t* sx_value, const Real_t* sy_value, const Real_t* sz_value, const Real_t* pole_value){

  size_t k0_max=1+sph_order;
  size_t k1_max=2*sph_order;
  size_t M_ves=k0_max*k1_max;

  Real_t dR2, dR2_p0, dR2_p1;
  { // dR2
    int k=k1_+k0_*k1_max;
    Real_t dR[COORD_DIM]={sx_coord[k]-tcoord[0],
                          sy_coord[k]-tcoord[1],
                          sz_coord[k]-tcoord[2]};
    dR2=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
  }
  { // dR2_p0
    Real_t dR[COORD_DIM]={pole_coord[0*COORD_DIM+0]-tcoord[0],
                          pole_coord[0*COORD_DIM+1]-tcoord[1],
                          pole_coord[0*COORD_DIM+2]-tcoord[2]};
    dR2_p0=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
  }
  { // dR2_p1
    Real_t dR[COORD_DIM]={pole_coord[1*COORD_DIM+0]-tcoord[0],
                          pole_coord[1*COORD_DIM+1]-tcoord[1],
                          pole_coord[1*COORD_DIM+2]-tcoord[2]};
    dR2_p1=dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2];
  }

  static pvfmm::Matrix<Real_t> M_interp;
  if(M_interp.Dim(0)!=k1_max || M_interp.Dim(1)!=8){ // Set M_interp
    #pragma omp critical
    if(M_interp.Dim(0)!=k1_max || M_interp.Dim(1)!=8){
      M_interp.ReInit(k1_max,8);
      M_interp.SetZero();
      size_t sph_order=k1_max/2;
      assert(2*sph_order==k1_max);
      for(size_t k0=0;k0<8;k0++)
      for(size_t k1=0;k1<k1_max;k1++)
      for(size_t k2=0;k2<=sph_order;k2++){
        M_interp[k1][k0]+=(!k2 || k2==sph_order?1.0:2.0)*sin(2.0*M_PI*k2*(k0*1.0/8.0))*sin(2.0*M_PI*k2*(k1*1.0/k1_max))*1.0/k1_max;
        M_interp[k1][k0]+=(!k2 || k2==sph_order?1.0:2.0)*cos(2.0*M_PI*k2*(k0*1.0/8.0))*cos(2.0*M_PI*k2*(k1*1.0/k1_max))*1.0/k1_max;
      }
    }
  }

  if(dR2_p0<dR2){ // Create patch
    Real_t tmp[3*8];
    pvfmm::Matrix<Real_t> value(3,8,tmp,false);
    { // Set value
      int k=0*k1_max;
      pvfmm::Matrix<Real_t> X_in, X_out;

      X_in .ReInit(1,k1_max,(Real_t*)&sx_value[k],false);
      X_out.ReInit(1,     8,             value[0],false);
      X_out=X_in*M_interp;

      X_in .ReInit(1,k1_max,(Real_t*)&sy_value[k],false);
      X_out.ReInit(1,     8,             value[1],false);
      X_out=X_in*M_interp;

      X_in .ReInit(1,k1_max,(Real_t*)&sz_value[k],false);
      X_out.ReInit(1,     8,             value[2],false);
      X_out=X_in*M_interp;
    }

    { // (0,0)
      patch_value_[(0*3+0)*COORD_DIM+0]+=value[0][0];
      patch_value_[(0*3+0)*COORD_DIM+1]+=value[1][0];
      patch_value_[(0*3+0)*COORD_DIM+2]+=value[2][0];
    }
    { // (0,1)
      patch_value_[(0*3+1)*COORD_DIM+0]+=value[0][7];
      patch_value_[(0*3+1)*COORD_DIM+1]+=value[1][7];
      patch_value_[(0*3+1)*COORD_DIM+2]+=value[2][7];
    }
    { // (0,2)
      patch_value_[(0*3+2)*COORD_DIM+0]+=value[0][6];
      patch_value_[(0*3+2)*COORD_DIM+1]+=value[1][6];
      patch_value_[(0*3+2)*COORD_DIM+2]+=value[2][6];
    }

    { // (1,0)
      patch_value_[(1*3+0)*COORD_DIM+0]+=value[0][1];
      patch_value_[(1*3+0)*COORD_DIM+1]+=value[1][1];
      patch_value_[(1*3+0)*COORD_DIM+2]+=value[2][1];
    }
    { // (1,1)
      patch_value_[(1*3+1)*COORD_DIM+0]+=pole_value[0*COORD_DIM+0];
      patch_value_[(1*3+1)*COORD_DIM+1]+=pole_value[0*COORD_DIM+1];
      patch_value_[(1*3+1)*COORD_DIM+2]+=pole_value[0*COORD_DIM+2];
    }
    { // (1,2)
      patch_value_[(1*3+2)*COORD_DIM+0]+=value[0][5];
      patch_value_[(1*3+2)*COORD_DIM+1]+=value[1][5];
      patch_value_[(1*3+2)*COORD_DIM+2]+=value[2][5];
    }

    { // (2,0)
      patch_value_[(2*3+0)*COORD_DIM+0]+=value[0][2];
      patch_value_[(2*3+0)*COORD_DIM+1]+=value[1][2];
      patch_value_[(2*3+0)*COORD_DIM+2]+=value[2][2];
    }
    { // (2,1)
      patch_value_[(2*3+1)*COORD_DIM+0]+=value[0][3];
      patch_value_[(2*3+1)*COORD_DIM+1]+=value[1][3];
      patch_value_[(2*3+1)*COORD_DIM+2]+=value[2][3];
    }
    { // (2,2)
      patch_value_[(2*3+2)*COORD_DIM+0]+=value[0][4];
      patch_value_[(2*3+2)*COORD_DIM+1]+=value[1][4];
      patch_value_[(2*3+2)*COORD_DIM+2]+=value[2][4];
    }
  }else if(dR2_p1<dR2){
    Real_t tmp[3*8];
    pvfmm::Matrix<Real_t> value(3,8,tmp,false);
    { // Set value
      int k=(k0_max-1)*k1_max;
      pvfmm::Matrix<Real_t> X_in, X_out;

      X_in .ReInit(1,k1_max,(Real_t*)&sx_value[k],false);
      X_out.ReInit(1,     8,             value[0],false);
      X_out=X_in*M_interp;

      X_in .ReInit(1,k1_max,(Real_t*)&sy_value[k],false);
      X_out.ReInit(1,     8,             value[1],false);
      X_out=X_in*M_interp;

      X_in .ReInit(1,k1_max,(Real_t*)&sz_value[k],false);
      X_out.ReInit(1,     8,             value[2],false);
      X_out=X_in*M_interp;
    }

    { // (0,0)
      patch_value_[(0*3+0)*COORD_DIM+0]+=value[0][0];
      patch_value_[(0*3+0)*COORD_DIM+1]+=value[1][0];
      patch_value_[(0*3+0)*COORD_DIM+2]+=value[2][0];
    }
    { // (0,1)
      patch_value_[(0*3+1)*COORD_DIM+0]+=value[0][1];
      patch_value_[(0*3+1)*COORD_DIM+1]+=value[1][1];
      patch_value_[(0*3+1)*COORD_DIM+2]+=value[2][1];
    }
    { // (0,2)
      patch_value_[(0*3+2)*COORD_DIM+0]+=value[0][2];
      patch_value_[(0*3+2)*COORD_DIM+1]+=value[1][2];
      patch_value_[(0*3+2)*COORD_DIM+2]+=value[2][2];
    }

    { // (1,0)
      patch_value_[(1*3+0)*COORD_DIM+0]+=value[0][7];
      patch_value_[(1*3+0)*COORD_DIM+1]+=value[1][7];
      patch_value_[(1*3+0)*COORD_DIM+2]+=value[2][7];
    }
    { // (1,1)
      patch_value_[(1*3+1)*COORD_DIM+0]+=pole_value[1*COORD_DIM+0];
      patch_value_[(1*3+1)*COORD_DIM+1]+=pole_value[1*COORD_DIM+1];
      patch_value_[(1*3+1)*COORD_DIM+2]+=pole_value[1*COORD_DIM+2];
    }
    { // (1,2)
      patch_value_[(1*3+2)*COORD_DIM+0]+=value[0][3];
      patch_value_[(1*3+2)*COORD_DIM+1]+=value[1][3];
      patch_value_[(1*3+2)*COORD_DIM+2]+=value[2][3];
    }

    { // (2,0)
      patch_value_[(2*3+0)*COORD_DIM+0]+=value[0][6];
      patch_value_[(2*3+0)*COORD_DIM+1]+=value[1][6];
      patch_value_[(2*3+0)*COORD_DIM+2]+=value[2][6];
    }
    { // (2,1)
      patch_value_[(2*3+1)*COORD_DIM+0]+=value[0][5];
      patch_value_[(2*3+1)*COORD_DIM+1]+=value[1][5];
      patch_value_[(2*3+1)*COORD_DIM+2]+=value[2][5];
    }
    { // (2,2)
      patch_value_[(2*3+2)*COORD_DIM+0]+=value[0][4];
      patch_value_[(2*3+2)*COORD_DIM+1]+=value[1][4];
      patch_value_[(2*3+2)*COORD_DIM+2]+=value[2][4];
    }
  }else{
    { // (0,0)
      if(k0_>0){
        int k=((k1_+k1_max-1)%k1_max)+(k0_-1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(0*3+0)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(0*3+0)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(0*3+0)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(0*3+0)*COORD_DIM+0]+=pole_value[0*COORD_DIM+0];
        patch_value_[(0*3+0)*COORD_DIM+1]+=pole_value[0*COORD_DIM+1];
        patch_value_[(0*3+0)*COORD_DIM+2]+=pole_value[0*COORD_DIM+2];
      }
    }
    { // (0,1)
      if(k0_>0){
        int k=  k1_                  +(k0_-1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(0*3+1)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(0*3+1)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(0*3+1)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(0*3+1)*COORD_DIM+0]+=pole_value[0*COORD_DIM+0];
        patch_value_[(0*3+1)*COORD_DIM+1]+=pole_value[0*COORD_DIM+1];
        patch_value_[(0*3+1)*COORD_DIM+2]+=pole_value[0*COORD_DIM+2];
      }
    }
    { // (0,2)
      if(k0_>0){
        int k=((k1_       +1)%k1_max)+(k0_-1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(0*3+2)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(0*3+2)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(0*3+2)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(0*3+2)*COORD_DIM+0]+=pole_value[0*COORD_DIM+0];
        patch_value_[(0*3+2)*COORD_DIM+1]+=pole_value[0*COORD_DIM+1];
        patch_value_[(0*3+2)*COORD_DIM+2]+=pole_value[0*COORD_DIM+2];
      }
    }

    { // (1,0)
      int k=((k1_+k1_max-1)%k1_max)+k0_*k1_max;
      assert(k>=0 && k<M_ves);
      patch_value_[(1*3+0)*COORD_DIM+0]+=sx_value[k];
      patch_value_[(1*3+0)*COORD_DIM+1]+=sy_value[k];
      patch_value_[(1*3+0)*COORD_DIM+2]+=sz_value[k];
    }
    { // (1,1)
      int k=  k1_                  +k0_*k1_max;
      assert(k>=0 && k<M_ves);
      patch_value_[(1*3+1)*COORD_DIM+0]+=sx_value[k];
      patch_value_[(1*3+1)*COORD_DIM+1]+=sy_value[k];
      patch_value_[(1*3+1)*COORD_DIM+2]+=sz_value[k];
    }
    { // (1,2)
      int k=((k1_       +1)%k1_max)+k0_*k1_max;
      assert(k>=0 && k<M_ves);
      patch_value_[(1*3+2)*COORD_DIM+0]+=sx_value[k];
      patch_value_[(1*3+2)*COORD_DIM+1]+=sy_value[k];
      patch_value_[(1*3+2)*COORD_DIM+2]+=sz_value[k];
    }

    { // (2,0)
      if(k0_<k0_max-1){
        int k=((k1_+k1_max-1)%k1_max)+(k0_+1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(2*3+0)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(2*3+0)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(2*3+0)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(2*3+0)*COORD_DIM+0]+=pole_value[1*COORD_DIM+0];
        patch_value_[(2*3+0)*COORD_DIM+1]+=pole_value[1*COORD_DIM+1];
        patch_value_[(2*3+0)*COORD_DIM+2]+=pole_value[1*COORD_DIM+2];
      }
    }
    { // (2,1)
      if(k0_<k0_max-1){
        int k=  k1_                  +(k0_+1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(2*3+1)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(2*3+1)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(2*3+1)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(2*3+1)*COORD_DIM+0]+=pole_value[1*COORD_DIM+0];
        patch_value_[(2*3+1)*COORD_DIM+1]+=pole_value[1*COORD_DIM+1];
        patch_value_[(2*3+1)*COORD_DIM+2]+=pole_value[1*COORD_DIM+2];
      }
    }
    { // (2,2)
      if(k0_<k0_max-1){
        int k=((k1_       +1)%k1_max)+(k0_+1)*k1_max;
        assert(k>=0 && k<M_ves);
        patch_value_[(2*3+2)*COORD_DIM+0]+=sx_value[k];
        patch_value_[(2*3+2)*COORD_DIM+1]+=sy_value[k];
        patch_value_[(2*3+2)*COORD_DIM+2]+=sz_value[k];
      }else{
        patch_value_[(2*3+2)*COORD_DIM+0]+=pole_value[1*COORD_DIM+0];
        patch_value_[(2*3+2)*COORD_DIM+1]+=pole_value[1*COORD_DIM+1];
        patch_value_[(2*3+2)*COORD_DIM+2]+=pole_value[1*COORD_DIM+2];
      }
    }
  }
}


template<typename Surf_t>
const NearSingular<Surf_t>::PVFMMVec_t& NearSingular<Surf_t>::operator()(bool update){
  if((update && update_interp) || (update_interp & NearSingular::UpdateSurfaceVel)){ // Compute near-singular
    pvfmm::Profile::Tic("NearInteraction",&comm,true);
    bool prof_state=pvfmm::Profile::Enable(false);
    size_t omp_p=omp_get_max_threads();
    SetupCoordData();

    Real_t&                  r_near=coord_setup.        r_near;
    PVFMMVec_t&           trg_coord=coord_setup.near_trg_coord;
    pvfmm::Vector<size_t>&  trg_cnt=coord_setup.  near_trg_cnt;
    pvfmm::Vector<size_t>&  trg_dsp=coord_setup.  near_trg_dsp;
    assert(S_vel);

    if(update)
    vel_interp.ReInit(trg_coord.Dim());
    vel_surfac.ReInit(trg_coord.Dim());

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

    pvfmm::Profile::Tic("Proj",&comm,true);
    size_t N_trg=trg_coord.Dim()/COORD_DIM;
    std::vector<char>   is_surf_pt      (N_trg);           // If a target point is a surface point
    std::vector<char>   is_extr_pt      (N_trg);           // If a target point is an exterior point
    std::vector<Real_t> proj_coord      (N_trg*COORD_DIM); // Projection coordinates (x,y,z)
    std::vector<Real_t> proj_veloc      (N_trg*COORD_DIM); // Velocity at projection coordinates
    std::vector<Real_t> proj_patch_param(N_trg*        2); // Projection patch prameter coordinates
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Setup
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
        // read each component of x
        const Real_t* sx_coord=x.getSubN_begin(i)+0*M_ves;
        const Real_t* sy_coord=x.getSubN_begin(i)+1*M_ves;
        const Real_t* sz_coord=x.getSubN_begin(i)+2*M_ves;

        Real_t pole_coord[2*COORD_DIM];
        { // Set pole values: pole_coord, pole_veloc
          for(size_t j=0;j<2*COORD_DIM;j++) pole_coord[j]=0;
          for(size_t k0=0;k0<k0_max;k0++){
            for(size_t k1=0;k1<k1_max;k1++){
              size_t k=k1+k0*k1_max;
              pole_coord[0*COORD_DIM+0]+=pole_quad[k0_max-1-k0]*sx_coord[k];
              pole_coord[0*COORD_DIM+1]+=pole_quad[k0_max-1-k0]*sy_coord[k];
              pole_coord[0*COORD_DIM+2]+=pole_quad[k0_max-1-k0]*sz_coord[k];
              pole_coord[1*COORD_DIM+0]+=pole_quad[         k0]*sx_coord[k];
              pole_coord[1*COORD_DIM+1]+=pole_quad[         k0]*sy_coord[k];
              pole_coord[1*COORD_DIM+2]+=pole_quad[         k0]*sz_coord[k];
            }
          }
        }

        // Determine if trg point is on the surface
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          size_t src_idx=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;

          is_surf_pt[trg_idx]=(trg_coord[trg_idx*COORD_DIM+0]==sx_coord[src_idx] &&
                               trg_coord[trg_idx*COORD_DIM+1]==sy_coord[src_idx] &&
                               trg_coord[trg_idx*COORD_DIM+2]==sz_coord[src_idx]);
        }

        // Compute projection for near points
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          if(is_surf_pt[trg_idx]) continue;

          QuadraticPatch<Real_t> patch_coord;
          { // Find nearest point on mesh and create patch
            int k0_, k1_; // mesh coordinates for nearest point
            { // Find nearest point on mesh
              size_t k=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
              k0_=k/k1_max;
              k1_=k%k1_max;
            }

            Real_t patch_coord_[3*3*COORD_DIM];
            { // compute patch_coord_
              for(size_t k=0;k<3*3*COORD_DIM;k++) patch_coord_[k]=0;
              patch_mesh<Real_t>(patch_coord_, k0_max-1, k0_, k1_,
                                 sx_coord, sy_coord, sz_coord, pole_coord, &trg_coord[trg_idx*COORD_DIM],
                                 sx_coord, sy_coord, sz_coord, pole_coord);
            }
            patch_coord=QuadraticPatch<Real_t>(&patch_coord_[0],COORD_DIM);
          }

          { // Find nearest point on patch (first interpolation point)
            Real_t& x=proj_patch_param[trg_idx*2+0];
            Real_t& y=proj_patch_param[trg_idx*2+1];
            is_extr_pt[trg_idx]=
            patch_coord.project(&  trg_coord[trg_idx*COORD_DIM],x,y);
            patch_coord.eval(x,y,&proj_coord[trg_idx*COORD_DIM]);
          }
        }


        if(!update || !update_interp) continue; ///////////////////////////////


        // read each component of S_vel
        const Real_t* sx_veloc=S_vel->getSubN_begin(i)+0*M_ves;
        const Real_t* sy_veloc=S_vel->getSubN_begin(i)+1*M_ves;
        const Real_t* sz_veloc=S_vel->getSubN_begin(i)+2*M_ves;

        // read each component of force_double
        const Real_t* sx_dforce=(force_double?force_double->getSubN_begin(i)+0*M_ves:NULL);
        const Real_t* sy_dforce=(force_double?force_double->getSubN_begin(i)+1*M_ves:NULL);
        const Real_t* sz_dforce=(force_double?force_double->getSubN_begin(i)+2*M_ves:NULL);

        Real_t pole_veloc[2*COORD_DIM];
        Real_t pole_force_dbl[2*COORD_DIM];
        { // Set pole values: pole_coord, pole_veloc
          for(size_t j=0;j<2*COORD_DIM;j++) pole_veloc[j]=0;
          for(size_t j=0;j<2*COORD_DIM;j++) pole_force_dbl[j]=0;
          for(size_t k0=0;k0<k0_max;k0++){
            for(size_t k1=0;k1<k1_max;k1++){
              size_t k=k1+k0*k1_max;
              pole_veloc[0*COORD_DIM+0]+=pole_quad[k0_max-1-k0]*sx_veloc[k];
              pole_veloc[0*COORD_DIM+1]+=pole_quad[k0_max-1-k0]*sy_veloc[k];
              pole_veloc[0*COORD_DIM+2]+=pole_quad[k0_max-1-k0]*sz_veloc[k];
              pole_veloc[1*COORD_DIM+0]+=pole_quad[         k0]*sx_veloc[k];
              pole_veloc[1*COORD_DIM+1]+=pole_quad[         k0]*sy_veloc[k];
              pole_veloc[1*COORD_DIM+2]+=pole_quad[         k0]*sz_veloc[k];

              if(force_double){
                pole_force_dbl[0*COORD_DIM+0]+=pole_quad[k0_max-1-k0]*sx_dforce[k];
                pole_force_dbl[0*COORD_DIM+1]+=pole_quad[k0_max-1-k0]*sy_dforce[k];
                pole_force_dbl[0*COORD_DIM+2]+=pole_quad[k0_max-1-k0]*sz_dforce[k];
                pole_force_dbl[1*COORD_DIM+0]+=pole_quad[         k0]*sx_dforce[k];
                pole_force_dbl[1*COORD_DIM+1]+=pole_quad[         k0]*sy_dforce[k];
                pole_force_dbl[1*COORD_DIM+2]+=pole_quad[         k0]*sz_dforce[k];
              }
            }
          }
        }

        // Compute projection for near points
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          if(is_surf_pt[trg_idx]) continue;

          QuadraticPatch<Real_t> patch_veloc;
          { // Find nearest point on mesh and create patch
            int k0_, k1_; // mesh coordinates for nearest point
            { // Find nearest point on mesh
              size_t k=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
              k0_=k/k1_max;
              k1_=k%k1_max;
            }

            Real_t patch_veloc_[3*3*COORD_DIM];
            { // compute patch_veloc_
              for(size_t k=0;k<3*3*COORD_DIM;k++) patch_veloc_[k]=0;
              patch_mesh<Real_t>(patch_veloc_, k0_max-1, k0_, k1_,
                                 sx_coord, sy_coord, sz_coord, pole_coord, &trg_coord[trg_idx*COORD_DIM],
                                 sx_veloc, sy_veloc, sz_veloc, pole_veloc);
            }
            if(force_double){ // add contribution from double layer to patch_veloc
              Real_t patch_force_dbl_[3*3*COORD_DIM];
              for(size_t k=0;k<3*3*COORD_DIM;k++) patch_force_dbl_[k]=0;
              patch_mesh<Real_t>(patch_force_dbl_, k0_max-1, k0_, k1_,
                                 sx_coord, sy_coord, sz_coord, pole_coord, &trg_coord[trg_idx*COORD_DIM],
                                 sx_dforce, sy_dforce, sz_dforce, pole_force_dbl);
              Real_t scal=0.5;
              if(is_extr_pt[trg_idx]) scal=-0.5;
              for(size_t k=0;k<3*3*COORD_DIM;k++) patch_veloc_[k]+=scal*patch_force_dbl_[k];
            }
            patch_veloc=QuadraticPatch<Real_t>(&patch_veloc_[0],COORD_DIM);
          }

          { // Find nearest point on patch (first interpolation point)
            Real_t& x=proj_patch_param[trg_idx*2+0];
            Real_t& y=proj_patch_param[trg_idx*2+1];
            patch_veloc.eval(x,y,&proj_veloc[trg_idx*COORD_DIM]);
          }
        }
      }
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("VelocInterp",&comm,true);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Compute vel_interp.
      PVFMMVec_t s_coord(M_ves*COORD_DIM);
      PVFMMVec_t interp_coord;
      PVFMMVec_t interp_veloc;
      PVFMMVec_t interp_x;

      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      if(update && update_interp)
      for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
        { // Set s_coord
          // read each component of x
          const Real_t* xk=x.getSubN_begin(i)+0*M_ves;
          const Real_t* yk=x.getSubN_begin(i)+1*M_ves;
          const Real_t* zk=x.getSubN_begin(i)+2*M_ves;
          for(size_t j=0;j<M_ves;j++){
            s_coord[j*COORD_DIM+0]=xk[j];
            s_coord[j*COORD_DIM+1]=yk[j];
            s_coord[j*COORD_DIM+2]=zk[j];
          }
        }

        { // Resize interp_coord, interp_veloc
          size_t near_cnt=0;
          for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
            size_t trg_idx=trg_dsp[i]+j;
            Real_t* t_veloc_j=&vel_interp[trg_idx*COORD_DIM];
            if(is_surf_pt[trg_idx]){
              t_veloc_j[0]=0.0;
              t_veloc_j[1]=0.0;
              t_veloc_j[2]=0.0;
            }else near_cnt++;
          }
          interp_coord.Resize(near_cnt*(INTERP_DEG-1)*COORD_DIM);
          interp_veloc.Resize(near_cnt*(INTERP_DEG-1)*COORD_DIM);
          interp_x.Resize(near_cnt);
          interp_veloc.SetZero();
        }

        { // Set interp_coord
          size_t near_cnt=0;
          for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
            size_t trg_idx=trg_dsp[i]+j;
            if(is_surf_pt[trg_idx]) continue;

            Real_t interp_coord0[COORD_DIM];
            interp_coord0[0]=proj_coord[trg_idx*COORD_DIM+0];
            interp_coord0[1]=proj_coord[trg_idx*COORD_DIM+1];
            interp_coord0[2]=proj_coord[trg_idx*COORD_DIM+2];

            { // Set interp_coord
              Real_t dR[COORD_DIM]={trg_coord[trg_idx*COORD_DIM+0]-interp_coord0[0],
                                    trg_coord[trg_idx*COORD_DIM+1]-interp_coord0[1],
                                    trg_coord[trg_idx*COORD_DIM+2]-interp_coord0[2]};
              Real_t dR_norm=sqrt(dR[0]*dR[0]+dR[1]*dR[1]+dR[2]*dR[2]);
              Real_t OOdR=1.0/dR_norm;

              interp_x[near_cnt]=dR_norm/r_near;
              for(size_t l=0;l<INTERP_DEG-1;l++){
                Real_t x=InterPoints<Real_t>(l+1, INTERP_DEG);
                interp_coord[(l+near_cnt*(INTERP_DEG-1))*COORD_DIM+0]=interp_coord0[0]+dR[0]*OOdR*r_near*x;
                interp_coord[(l+near_cnt*(INTERP_DEG-1))*COORD_DIM+1]=interp_coord0[1]+dR[1]*OOdR*r_near*x;
                interp_coord[(l+near_cnt*(INTERP_DEG-1))*COORD_DIM+2]=interp_coord0[2]+dR[2]*OOdR*r_near*x;
              }
            }
            near_cnt++;
          }
        }

        if(interp_veloc.Dim()){ // Compute velocity at interpolation points
          if(qforce_single){
            stokes_sl(&s_coord[0], M_ves, &qforce_single[0][0]+M_ves*(COORD_DIM*1)*i, 1, &interp_coord[0], interp_coord.Dim()/COORD_DIM, &interp_veloc[0], NULL);
          }
          if(qforce_double){
            stokes_dl(&s_coord[0], M_ves, &qforce_double[0][0]+M_ves*(COORD_DIM*2)*i, 1, &interp_coord[0], interp_coord.Dim()/COORD_DIM, &interp_veloc[0], NULL);
          }
        }

        { // Interpolate
          size_t near_cnt=0;
          for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
            size_t trg_idx=trg_dsp[i]+j;
            if(is_surf_pt[trg_idx]) continue;
            Real_t* t_veloc_j=&vel_interp[trg_idx*COORD_DIM];

            { // Interpolate and find target velocity t_veloc_j
              pvfmm::Matrix<Real_t> y(INTERP_DEG,1);
              pvfmm::Matrix<Real_t> x(INTERP_DEG,1);
              Real_t x_=interp_x[near_cnt];
              for(size_t l=0;l<INTERP_DEG;l++){
                x[l][0]=InterPoly(x_,l,INTERP_DEG);
              }
              for(size_t k=0;k<COORD_DIM;k++){
                y[0][0]=proj_veloc[trg_idx*COORD_DIM+k];
                for(size_t l=0;l<INTERP_DEG-1;l++){
                  y[l+1][0]=interp_veloc[(l+near_cnt*(INTERP_DEG-1))*COORD_DIM+k];
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
                t_veloc_j[k]=0;
                for(size_t l=0;l<INTERP_DEG;l++) t_veloc_j[k]+=y[l][0]*x[l][0];
              }
            }
            near_cnt++;
          }
        }
      }
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("VelocSurf",&comm,true);
    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p;tid++){ // Compute vel_surfac.
      size_t a=((tid+0)*N_ves)/omp_p;
      size_t b=((tid+1)*N_ves)/omp_p;
      for(size_t i=a;i<b;i++) if(trg_cnt[i]){ // loop over all vesicles
        // read each component of veloc
        const Real_t* sx_veloc=S_vel->getSubN_begin(i)+0*M_ves;
        const Real_t* sy_veloc=S_vel->getSubN_begin(i)+1*M_ves;
        const Real_t* sz_veloc=S_vel->getSubN_begin(i)+2*M_ves;

        // Get velocity for surface points
        for(size_t j=0;j<trg_cnt[i];j++){ // loop over target points
          size_t trg_idx=trg_dsp[i]+j;
          if(is_surf_pt[trg_idx]){
            size_t src_idx=coord_setup.near_ves_pt_id[trg_idx]-M_ves*i;
            vel_surfac[trg_idx*COORD_DIM+0]=sx_veloc[src_idx];
            vel_surfac[trg_idx*COORD_DIM+1]=sy_veloc[src_idx];
            vel_surfac[trg_idx*COORD_DIM+2]=sz_veloc[src_idx];
          }else{
            vel_surfac[trg_idx*COORD_DIM+0]=0.0;
            vel_surfac[trg_idx*COORD_DIM+1]=0.0;
            vel_surfac[trg_idx*COORD_DIM+2]=0.0;
          }
        }
      }
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Tic("Scatter",&comm,true);
    if(update)
    VelocityScatter(vel_interp);
    VelocityScatter(vel_surfac);

    assert(vel_interp.Dim()==T.Dim());
    assert(vel_surfac.Dim()==T.Dim());
    #pragma omp parallel for
    for(size_t i=0;i<vel_surfac.Dim();i++){
      vel_surfac[i]+=vel_interp[i];
    }

    if(update){
      update_interp=NearSingular::UpdateNone;
    }else{
      update_interp=(update_interp & ~(NearSingular::UpdateSurfaceVel));
    }
    pvfmm::Profile::Toc();

    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
  }
  return vel_surfac;
}

#endif
