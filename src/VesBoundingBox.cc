#include <omp.h>
#include <iostream>

#include <ompUtils.h>
#include <parUtils.h>
#include <profile.hpp>
#include <mortonid.hpp>
#include <SphericalHarmonics.h>
#include <mpi_tree.hpp> // Only for vis

template<typename Real_t>
VesBoundingBox<Real_t>::VesBoundingBox(Real_t box_size, MPI_Comm c){
    box_size_=box_size;
    comm=c;
}

template<typename Real_t>
void VesBoundingBox<Real_t>::SetVesBoundingBox(const PVFMMVec_t& ves_coord_s, const PVFMMVec_t& ves_coord_e,
        const int ves_stride, Real_t min_sep)
{
    // calculate the bounding boxes for vesicles with start, end position and min_sep
}

template<typename Real_t>
template<typename Vec>
void VesBoundingBox<Real_t>::SetVesBoundingBox(const Vec& ves_coord_s, const Vec& ves_coord_e, Real_t min_sep)
{
    // calculate the bounding boxes for vesicles with start, end position and min_sep
}

template<typename Real_t>
void VesBoundingBox<Real_t>::SetVesBoundingBox(const PVFMMVec_t& BB_min, const PVFMMVec_t& BB_max)
{
    BB_min_ = BB_min;
    BB_max_ = BB_max;
}

template<typename Real_t>
void VesBoundingBox<Real_t>::GetContactBoundingBoxPair(std::vector< std::pair<size_t, size_t> > &BBIPairs)
{
    assert(BB_pts_); assert(BB_min_); assert(BB_max_);
    
    MPI_Comm_size(comm, &np_);
    MPI_Comm_rank(comm, &rank_);
    omp_p_ = omp_get_max_threads();

#ifdef VES_BOUNDING_BOX_DEBUG
    INFO("rank: "<<rank_<<" of "<<np_<<" processes has "<<omp_p_<<" threads");
#endif

    N_bbox_ = BB_pts_.Dim()/COORD_DIM;
    assert(N_bbox_ > 0);
    assert(N_bbox_*COORD_DIM == BB_pts_.Dim()); 
    assert(BB_pts_.Dim() == BB_min_.Dim()); assert(BB_pts_.Dim() == BB_max_.Dim());

    pvfmm::Profile::Tic("GetContactBBPair", &comm, true);
    bool prof_state=pvfmm::Profile::Enable(false);
    
    TREEGRID BB_let;

    {   // construct BB_let
        pvfmm::Profile::Tic("BBLET", &comm, true);

        SetTreeParams();
#ifdef VES_BOUNDING_BOX_DEBUG
        INFO("bbox_: "<<bbox_[0]<<","<<bbox_[1]<<","<<bbox_[2]<<","<<bbox_[3]<<".");
        INFO("r_near_: "<<r_near_);
        INFO("tree_depth_: "<<tree_depth_);
#endif
        
        GenerateBBPoints();
        
        ConstructLocalTree(BB_let);
#ifdef VES_BOUNDING_BOX_DEBUG
        INFO("size of mid: "<<BB_let.mid.Dim());
        INFO(BB_let.mid);
        INFO("size of pt_cnt: "<<BB_let.pt_cnt.Dim());
        INFO(BB_let.pt_cnt);
        INFO("size of pt_dsp: "<<BB_let.pt_dsp.Dim());
        INFO(BB_let.pt_dsp);
        INFO("size of mins: "<<BB_let.mins.Dim());
        INFO(BB_let.mins);
        //INFO("size of pt_coord: "<<BB_let.pt_coord.Dim());
        //INFO(BB_let.pt_coord);
        INFO("size of pt_id: "<<BB_let.pt_id.Dim());
        INFO(BB_let.pt_id);
#endif
        
        pvfmm::Profile::Toc();
    }

    {   // find contact BB pair
        pvfmm::Vector<size_t> near_trg_box_id;
        pvfmm::Vector<size_t> near_src_box_id;
        FindNearPair(BB_let, near_trg_box_id, near_src_box_id);
#ifdef VES_BOUNDING_BOX_DEBUG
        INFO("size of near_trg_box_id: "<<near_trg_box_id.Dim());
        INFO(near_trg_box_id);
        INFO("size of near_src_pt_id: "<<near_src_box_id.Dim());
        INFO(near_src_box_id);
#endif
    }

    pvfmm::Profile::Enable(prof_state);
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void VesBoundingBox<Real_t>::SetTreeParams()
{
    pvfmm::Profile::Tic("TreeParams", &comm, true);

    Real_t* bbox = bbox_;
    Real_t& r_near = r_near_;
    size_t& tree_depth = tree_depth_;

    tree_depth = 0;
    Real_t r_bbox_max = 0;
    Real_t r_bbox_min = std::numeric_limits<Real_t>::max();
    Real_t r_bbox_avg = 0;
    assert(N_bbox_ > 0);
    {   // determine r_bbox
        std::vector<Real_t> r2_bbox_max_omp(omp_p_);
        std::vector<Real_t> r2_bbox_min_omp(omp_p_);
        std::vector<Real_t> r2_bbox_avg_omp(omp_p_);
        #pragma omp parallel for
        for(size_t tid=0; tid<omp_p_; tid++){
            size_t a = ((tid+0)*N_bbox_)/omp_p_;
            size_t b = ((tid+1)*N_bbox_)/omp_p_;

            Real_t r2_bbox_max = 0;
            Real_t r2_bbox_min = std::numeric_limits<Real_t>::max();
            Real_t r2_bbox_avg = 0;
            for(size_t i=a; i<b; i++){ // compute r2_bbox
                const Real_t* min_i = &BB_min_[i*COORD_DIM];
                const Real_t* max_i = &BB_max_[i*COORD_DIM];

                Real_t dx = min_i[0] - max_i[0];
                Real_t dy = min_i[1] - max_i[1];
                Real_t dz = min_i[2] - max_i[2];
                Real_t r2 = dx*dx + dy*dy + dz*dz;

                r2_bbox_max = std::max(r2_bbox_max, r2);
                r2_bbox_min = std::min(r2_bbox_min, r2);
            }
            r2_bbox_max_omp[tid] = r2_bbox_max;
            r2_bbox_min_omp[tid] = r2_bbox_min;
        }

        {   // determine r_near (global max)
            double r_bbox_max_loc = 0, r_bbox_max_glb = 0;
            double r_bbox_min_loc = 0, r_bbox_min_glb = 0;
            for(size_t tid=0;tid<omp_p_;tid++){
                r_bbox_max_loc = std::max(r2_bbox_max_omp[tid], r_bbox_max_loc);
                r_bbox_min_loc = std::min(r2_bbox_min_omp[tid], r_bbox_min_loc);
            }
            r_bbox_max_loc = sqrt(r_bbox_max_loc);
            r_bbox_min_loc = sqrt(r_bbox_min_loc);
            MPI_Allreduce(&r_bbox_max_loc, &r_bbox_max_glb, 1, MPI_DOUBLE, MPI_MAX, comm);
            MPI_Allreduce(&r_bbox_min_loc, &r_bbox_min_glb, 1, MPI_DOUBLE, MPI_MIN, comm);
            r_bbox_max = r_bbox_max_glb;
            r_bbox_min = r_bbox_min_glb;
        }
    }

    r_near = r_bbox_min;
    if(box_size_>0 && r_bbox_max > box_size_){
        COUTDEBUG("Domain too small for bounding box size");
        assert(false);
        exit(0);
    }

    if(box_size_<=0){ // determine the global bbox, tree_depth
        Real_t scale_x, shift_x[COORD_DIM];
        Real_t scale_tmp;
        
        // determine bounding box
        GlobalBoundingBox(&scale_tmp, shift_x);

        { // scale_x, pt_tree_depth
            assert(scale_tmp!=0);
            
            Real_t domain_length = 1.0/scale_tmp + 4*r_near;
            INFO("domain_length_tmp: "<<domain_length);
            Real_t leaf_size = r_near;
            scale_x = 1.0/leaf_size;
            while(domain_length*scale_x>1.0 && tree_depth<MAX_DEPTH-1){
                scale_x *= 0.5;
                tree_depth++;
            }
            INFO("domain_length_actual: "<<1.0/scale_x);
            leaf_size_ = leaf_size;
        }

        for(size_t j=0;j<COORD_DIM;j++){ // Update shift_x
            shift_x[j]=((shift_x[j]/scale_tmp)+2*r_near)*scale_x;
        }

        bbox[0]=shift_x[0];
        bbox[1]=shift_x[1];
        bbox[2]=shift_x[2];
        bbox[3]=scale_x;
    }else{
        bbox[0] = 0;
        bbox[1] = 0;
        bbox[2] = 0;
        bbox[3] = 1.0/box_size_;

        // determine the tree depth
        Real_t leaf_size = box_size_;
        // 1*r_near < leaf_size <= 2*r_near
        while(leaf_size*0.5>r_near && tree_depth<MAX_DEPTH-1){
            leaf_size *= 0.5;
            tree_depth++;
        }
        leaf_size *= 0.5; tree_depth++;
        leaf_size_ = leaf_size;
    }
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void VesBoundingBox<Real_t>::ConstructLocalTree(TREEGRID &BB_let)
{
    pvfmm::Profile::Tic("LocalTree",&comm,true);
    //PVFMMVec_t           & pt_coord=BB_let.pt_coord;
    PVFMMVec_t           & box_min =BB_let.box_min;
    PVFMMVec_t           & box_max =BB_let.box_max;
    pvfmm::Vector<size_t>& pt_id   =BB_let.pt_id   ;
    pvfmm::Vector<size_t>& box_id  =BB_let.box_id   ;

    pvfmm::Vector<pvfmm::MortonId>& let_mid   =BB_let.mid;
    pvfmm::Vector<size_t>&          let_pt_cnt=BB_let.pt_cnt;
    pvfmm::Vector<size_t>&          let_pt_dsp=BB_let.pt_dsp;
    pvfmm::Vector<pvfmm::MortonId>& let_mins  =BB_let.mins;

    { // build scatter-indices (pt_id) and tree (let_mid, let_pt_cnt, let_pt_dsp)
        pvfmm::Vector<pvfmm::MortonId> pt_mid(N_bbox_);
        { // build pt_mid
            Real_t scale_x, shift_x[COORD_DIM];
            { // set scale_x, shift_x
                shift_x[0]=bbox_[0];
                shift_x[1]=bbox_[1];
                shift_x[2]=bbox_[2];
                scale_x=bbox_[3];
            }

            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                Real_t c[COORD_DIM];
                size_t a=((tid+0)*N_bbox_)/omp_p_;
                size_t b=((tid+1)*N_bbox_)/omp_p_;
                for(size_t i=a;i<b;i++){
                    for(size_t k=0;k<COORD_DIM;k++){
                        c[k]=BB_pts_[i*COORD_DIM+k]*scale_x+shift_x[k];
                        while(c[k]< 0.0) c[k]+=1.0;
                        while(c[k]>=1.0) c[k]-=1.0;
                    }
                    pt_mid[i]=pvfmm::MortonId(c,tree_depth_);
                }
            }
        }

        pt_id   .ReInit(N_bbox_);
        pvfmm::par::SortScatterIndex(pt_mid, pt_id, comm);
        pvfmm::par::ScatterForward  (pt_mid, pt_id, comm);
        { // build let_mins
            let_mins.ReInit(np_);
            MPI_Allgather(&  pt_mid[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                          &let_mins[0], 1, pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
            if(rank_) assert(let_mins[rank_]!=let_mins[rank_-1]);
            let_mins[0]=pvfmm::MortonId(0,0,0,tree_depth_);
        }
        { // Exchange shared octant with neighbour
            int send_size=0;
            int recv_size=0;
            if(rank_<np_-1){ // send_size
                send_size=pt_mid.Dim()-(std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), let_mins[rank_+1])-&pt_mid[0]);
            }
            { // recv_size
                MPI_Status status;
                MPI_Sendrecv(&send_size,1,MPI_INT,(rank_<np_-1?rank_+1:0),0,
                             &recv_size,1,MPI_INT,(rank_>0?rank_-1:np_-1),0,comm,&status);
            }

            { // Set new pt_id
                pvfmm::Vector<size_t> pt_id_new(pt_id.Dim()+recv_size-send_size);
                memcpy(&pt_id_new[0]+recv_size, &pt_id[0], (pt_id.Dim()-send_size)*sizeof(size_t));

                MPI_Status status;
                MPI_Sendrecv(&pt_id[0]+pt_id.Dim()-send_size,send_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank_<np_-1?rank_+1:0),0,
                             &pt_id_new[0]                  ,recv_size,pvfmm::par::Mpi_datatype<size_t>::value(),(rank_>0?rank_-1:np_-1),0,comm,&status);
                pt_id.Swap(pt_id_new);
            }
            { // Set new pt_mid
                pvfmm::Vector<pvfmm::MortonId> pt_mid_new(pt_mid.Dim()+recv_size-send_size);
                memcpy(&pt_mid_new[0]+recv_size, &pt_mid[0], (pt_mid.Dim()-send_size)*sizeof(pvfmm::MortonId));
                for(size_t i=0;i<recv_size;i++) pt_mid_new[i]=let_mins[rank_];
                pt_mid.Swap(pt_mid_new);
            }
        }
        { // Sort points by pt_id in each octant
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                size_t a=(pt_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(pt_mid.Dim()*(tid+1))/omp_p_;
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
            std::vector<std::vector<pvfmm::MortonId> > mid_omp(omp_p_);
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                std::vector<pvfmm::MortonId>& mid=mid_omp[tid];
                size_t a=(pt_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(pt_mid.Dim()*(tid+1))/omp_p_;
                if(a>0           ) a=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[a])-&pt_mid[0];
                if(b<pt_mid.Dim()) b=std::lower_bound(&pt_mid[0], &pt_mid[0]+pt_mid.Dim(), pt_mid[b])-&pt_mid[0];
                if(a<b) mid.push_back(pt_mid[a]);
                for(size_t i=a;i<b;i++){
                    if(mid.back()!=pt_mid[i]) mid.push_back(pt_mid[i]);
                }
            }
            { // Resize let_mid
                size_t size=0;
                for(size_t tid=0;tid<omp_p_;tid++){
                    size+=mid_omp[tid].size();
                }
                let_mid.ReInit(size);
            }
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){ // Set let_mid
                size_t offset=0;
                for(size_t i=0;i<tid;i++){
                    offset+=mid_omp[i].size();
                }

                std::vector<pvfmm::MortonId>& mid=mid_omp[tid];
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
    { // scatter pt_coord, box_id, box_min, box_max
        //pt_coord=BB_pts_;
        //pvfmm::par::ScatterForward(pt_coord, pt_id, comm);
        box_min=BB_pts_min_;
        pvfmm::par::ScatterForward(box_min, pt_id, comm);
        box_max=BB_pts_max_;
        pvfmm::par::ScatterForward(box_max, pt_id, comm);
        box_id = BB_id_;
        pvfmm::par::ScatterForward(box_id, pt_id, comm);
    }
    pvfmm::Profile::Toc();
}

/*
template<typename Real_t>
void VesBoundingBox<Real_t>::AddGhostNodes(TREEGRID &BB_let)
{
    pvfmm::Profile::Tic("AddGhostNodes");
    // TODO: Replace MPI_Alltoallv with Mpi_Alltoallv_sparse

    pvfmm::Vector<size_t> snode_id;
    pvfmm::Vector<int> snode_cnt(np_);
    pvfmm::Vector<int> snode_dsp(np_);
    pvfmm::Vector<int> rnode_cnt(np_);
    pvfmm::Vector<int> rnode_dsp(np_);
    {   // compute shared_pair
        pvfmm::Vector<pvfmm::par::SortPair<int,size_t> > shared_pair; // pid, node_id list
        { // Compute shared_pair
            std::vector<std::vector<pvfmm::par::SortPair<int,size_t> > > shared_pair_omp(omp_p_);
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                pvfmm::Vector<pvfmm::MortonId>&  let_mid=BB_let. mid;
                pvfmm::Vector<pvfmm::MortonId>& let_mins=BB_let.mins;

                Real_t coord[COORD_DIM];
                Real_t s=pow(0.5,tree_depth_);
                pvfmm::par::SortPair<int, size_t> pair;
                size_t a=(let_mid.Dim()*(tid+0))/omp_p_;
                size_t b=(let_mid.Dim()*(tid+1))/omp_p_;
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
                                        pvfmm::MortonId m(c,tree_depth_+1);
                                        if((rank_ && m<let_mins[rank_]) || (rank_<np_-1 && m>=let_mins[rank_+1])){
                                            int pid=std::lower_bound(&let_mins[0], &let_mins[0]+np_, m)-&let_mins[0]-1;
                                            assert(pid!=rank_);
                                            assert(pid>=0);
                                            assert(pid<np_);
                                            pid_set.insert(pid);
                                        }
                                    }
                                }
                    for(std::set<int>::iterator it=pid_set.begin(); it!=pid_set.end(); ++it){ // Add shared pair
                        pair.data=i;
                        pair.key=*it;
                        shared_pair_omp[tid].push_back(pair);
                    }
                }
            }
            { // Resize shared_pair
                size_t size=0;
                for(size_t tid=0;tid<omp_p_;tid++){
                    size+=shared_pair_omp[tid].size();
                }
                shared_pair.ReInit(size);
            }
            #pragma omp parallel for
            for(size_t tid=0;tid<omp_p_;tid++){
                size_t offset=0;
                for(size_t i=0;i<tid;i++){
                    offset+=shared_pair_omp[i].size();
                }
                for(size_t i=0;i<shared_pair_omp[tid].size();i++){
                    shared_pair[offset+i]=shared_pair_omp[tid][i];
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
            for(size_t tid=0;tid<omp_p_;tid++){
                size_t a=(tid+0)*shared_pair.Dim()/omp_p_;
                size_t b=(tid+1)*shared_pair.Dim()/omp_p_;
                size_t pid_a=(a<shared_pair.Dim()?shared_pair[a].key:np_);
                size_t pid_b=(b<shared_pair.Dim()?shared_pair[b].key:np_);
                for(size_t i=pid_a;i<pid_b;i++){
                    a=(i+0<np_?snode_dsp[i+0]:shared_pair.Dim());
                    b=(i+1<np_?snode_dsp[i+1]:shared_pair.Dim());
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
    {   // send-recv node data
        size_t send_size=snode_cnt[np_-1]+snode_dsp[np_-1];
        send_mid.ReInit(send_size);
        send_pt_cnt.ReInit(send_size);
        send_pt_dsp.ReInit(send_size+1);

        size_t recv_size=rnode_cnt[np_-1]+rnode_dsp[np_-1];
        recv_mid.ReInit(recv_size);
        recv_pt_cnt.ReInit(recv_size);
        recv_pt_dsp.ReInit(recv_size+1);

        for(size_t i=0;i<send_size;i++){ // Set send data
            send_mid   [i]=BB_let.mid   [snode_id[i]];
            send_pt_cnt[i]=BB_let.pt_cnt[snode_id[i]];
        }
        MPI_Alltoallv(&send_mid[0], &snode_cnt[0], &snode_dsp[0], pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(),
                      &recv_mid[0], &rnode_cnt[0], &rnode_dsp[0], pvfmm::par::Mpi_datatype<pvfmm::MortonId>::value(), comm);
        MPI_Alltoallv(&send_pt_cnt[0], &snode_cnt[0], &snode_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(),
                      &recv_pt_cnt[0], &rnode_cnt[0], &rnode_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(), comm);

        send_pt_dsp[0]=0; pvfmm::omp_par::scan(&send_pt_cnt[0], &send_pt_dsp[0], send_pt_cnt.Dim()+1);
        recv_pt_dsp[0]=0; pvfmm::omp_par::scan(&recv_pt_cnt[0], &recv_pt_dsp[0], recv_pt_cnt.Dim()+1);
    }

    PVFMMVec_t recv_pt_coord;
    pvfmm::Vector<size_t> recv_pt_id;

    PVFMMVec_t send_pt_coord;
    pvfmm::Vector<size_t> send_pt_id;
    {   // send-recv pt data
        size_t send_size_pt=send_pt_dsp[send_mid.Dim()];
        send_pt_coord.ReInit(send_size_pt*COORD_DIM);
        send_pt_id   .ReInit(send_size_pt);
        { // Set send data
            #pragma omp parallel for
            for(size_t i=0;i<send_pt_cnt.Dim();i++){
                size_t offset_in=BB_let.pt_dsp[snode_id[i]];
                size_t offset_out=send_pt_dsp[i];
                for(size_t j=0;j<send_pt_cnt[i];j++){
                    send_pt_coord[(offset_out+j)*COORD_DIM+0]=BB_let.pt_coord[(offset_in+j)*COORD_DIM+0];
                    send_pt_coord[(offset_out+j)*COORD_DIM+1]=BB_let.pt_coord[(offset_in+j)*COORD_DIM+1];
                    send_pt_coord[(offset_out+j)*COORD_DIM+2]=BB_let.pt_coord[(offset_in+j)*COORD_DIM+2];
                    send_pt_id   [ offset_out+j             ]=BB_let.pt_id   [ offset_in+j             ];
                }
            }
        }

        size_t recv_size_pt=recv_pt_dsp[recv_mid.Dim()];
        recv_pt_coord.ReInit(recv_size_pt*COORD_DIM);
        recv_pt_id   .ReInit(recv_size_pt);

        { // Send-recv data
            pvfmm::Vector<int> send_cnt(np_), send_dsp(np_);
            pvfmm::Vector<int> recv_cnt(np_), recv_dsp(np_);
            { // Set send_dsp, send_cnt,  recv_dsp, recv_cnt
                #pragma omp parallel for
                for(size_t i=0;i<np_;i++){
                    send_dsp[i]=send_pt_dsp[snode_dsp[i]];
                    recv_dsp[i]=recv_pt_dsp[rnode_dsp[i]];
                }
                #pragma omp parallel for
                for(size_t i=0;i<np_-1;i++){
                    send_cnt[i]=send_dsp[i+1]-send_dsp[i];
                    recv_cnt[i]=recv_dsp[i+1]-recv_dsp[i];
                }
                send_cnt[np_-1]=send_size_pt-send_dsp[np_-1];
                recv_cnt[np_-1]=recv_size_pt-recv_dsp[np_-1];
            }

            MPI_Alltoallv(&send_pt_id   [0], &send_cnt[0], &send_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(),
                          &recv_pt_id   [0], &recv_cnt[0], &recv_dsp[0], pvfmm::par::Mpi_datatype<size_t>::value(), comm);
            for(size_t i=0;i<np_;i++){
                send_cnt[i]*=COORD_DIM; send_dsp[i]*=COORD_DIM;
                recv_cnt[i]*=COORD_DIM; recv_dsp[i]*=COORD_DIM;
            }
            MPI_Alltoallv(&send_pt_coord[0], &send_cnt[0], &send_dsp[0], pvfmm::par::Mpi_datatype<Real_t>::value(),
                          &recv_pt_coord[0], &recv_cnt[0], &recv_dsp[0], pvfmm::par::Mpi_datatype<Real_t>::value(), comm);
        }
    }

    {   // add ghost nodes to BB_let
        pvfmm::Vector<pvfmm::MortonId> new_mid(BB_let.mid .Dim()+recv_mid     .Dim());
        pvfmm::Vector<size_t> new_pt_cnt  (BB_let.pt_cnt  .Dim()+recv_pt_cnt  .Dim());
        pvfmm::Vector<size_t> new_pt_dsp  (BB_let.pt_dsp  .Dim()+recv_pt_dsp  .Dim());

        PVFMMVec_t new_pt_coord(BB_let.pt_coord.Dim()+recv_pt_coord.Dim());
        pvfmm::Vector<size_t> new_pt_id   (BB_let.pt_id   .Dim()+recv_pt_id   .Dim());

        { // Copy mid
            size_t offset=0, size=0;
            size_t recv_split=rnode_dsp[rank_];
            size=                recv_split; memcpy(&new_mid[0]+offset, & recv_mid[0]           , size*sizeof(pvfmm::MortonId)); offset+=size;
            size=BB_let.mid.Dim()           ; memcpy(&new_mid[0]+offset, &BB_let.mid[0]           , size*sizeof(pvfmm::MortonId)); offset+=size;
            size= recv_mid.Dim()-recv_split; memcpy(&new_mid[0]+offset, & recv_mid[0]+recv_split, size*sizeof(pvfmm::MortonId));
        }
        { // Copy pt_cnt
            size_t offset=0, size=0;
            size_t recv_split=rnode_dsp[rank_];
            size=                   recv_split; memcpy(&new_pt_cnt[0]+offset, & recv_pt_cnt[0]           , size*sizeof(size_t)); offset+=size;
            size=BB_let.pt_cnt.Dim()           ; memcpy(&new_pt_cnt[0]+offset, &BB_let.pt_cnt[0]           , size*sizeof(size_t)); offset+=size;
            size= recv_pt_cnt.Dim()-recv_split; memcpy(&new_pt_cnt[0]+offset, & recv_pt_cnt[0]+recv_split, size*sizeof(size_t));
        }
        { // Compute pt_dsp
            new_pt_dsp[0]=0; pvfmm::omp_par::scan(&new_pt_cnt[0], &new_pt_dsp[0], new_pt_cnt.Dim());
        }

        { // Copy pt_coord
            size_t offset=0, size=0;
            size_t recv_split=recv_pt_dsp[rnode_dsp[rank_]]*COORD_DIM;
            size=                     recv_split; memcpy(&new_pt_coord[0]+offset, & recv_pt_coord[0]           , size*sizeof(Real_t)); offset+=size;
            size=BB_let.pt_coord.Dim()           ; memcpy(&new_pt_coord[0]+offset, &BB_let.pt_coord[0]           , size*sizeof(Real_t)); offset+=size;
            size= recv_pt_coord.Dim()-recv_split; memcpy(&new_pt_coord[0]+offset, & recv_pt_coord[0]+recv_split, size*sizeof(Real_t));
        }
        { // Copy pt_id
            size_t offset=0, size=0;
            size_t recv_split=recv_pt_dsp[rnode_dsp[rank_]];
            size=                  recv_split; memcpy(&new_pt_id[0]+offset, & recv_pt_id[0]           , size*sizeof(size_t)); offset+=size;
            size=BB_let.pt_id.Dim()           ; memcpy(&new_pt_id[0]+offset, &BB_let.pt_id[0]           , size*sizeof(size_t)); offset+=size;
            size= recv_pt_id.Dim()-recv_split; memcpy(&new_pt_id[0]+offset, & recv_pt_id[0]+recv_split, size*sizeof(size_t));
        }

        new_mid     .Swap(BB_let.mid     );
        new_pt_cnt  .Swap(BB_let.pt_cnt  );
        new_pt_dsp  .Swap(BB_let.pt_dsp  );

        new_pt_coord.Swap(BB_let.pt_coord);
        new_pt_id   .Swap(BB_let.pt_id   );
    }
    pvfmm::Profile::Toc();
}

template<typename Real_t>
void VesBoundingBox<Real_t>::FindNearPair(TREEGRID &BB_let, 
        pvfmm::Vector<size_t> &near_trg_pt_id, pvfmm::Vector<size_t> &near_src_pt_id)
{
    pvfmm::Profile::Tic("NearPair",&comm,true);
    PVFMMVec_t &pt_coord = BB_let.pt_coord;
    pvfmm::Vector<size_t> &pt_id = BB_let.pt_id;
    pvfmm::Vector<pvfmm::MortonId> &tree_mid = BB_let.mid;
    pvfmm::Vector<size_t> &tree_pt_cnt = BB_let.pt_cnt;
    pvfmm::Vector<size_t> &tree_pt_dsp = BB_let.pt_dsp;

    std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_omp(omp_p_); // (loc_trg_id, BB_let.pt_id)
    #pragma omp parallel num_threads(omp_p_)
    {
        size_t tid=omp_get_thread_num();
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];
        size_t tree_depth; Real_t r2_near;
        { // Set tree_depth, r_near
            tree_depth=BB_let.mid[0].GetDepth();
            r2_near=r_near_;
            r2_near*=r2_near;
        }
        Real_t s=pow(0.5,tree_depth);

        size_t FLOP=0;
        size_t a=((tid+0)*tree_mid.Dim())/omp_p_;
        size_t b=((tid+1)*tree_mid.Dim())/omp_p_;
        for(size_t i=a;i<b;i++){
            size_t tcnt=tree_pt_cnt[i];
            size_t tdsp=tree_pt_dsp[i];
            PVFMMVec_t tcoord;
            pvfmm::Vector<size_t> tpt_id;
            if(tcnt){ // Set t_coord
                tcoord.ReInit(tcnt*COORD_DIM,&pt_coord[tdsp*COORD_DIM],false);
                tpt_id.ReInit(tcnt,          &pt_id[tdsp]             ,false);
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
                                pvfmm::Vector<pvfmm::MortonId>& mid=BB_let.mid;
                                int k=std::lower_bound(&mid[0], &mid[0]+mid.Dim(), m)-&mid[0];
                                if(k<mid.Dim() && mid[k]==m){
                                    scoord[indx].ReInit(BB_let.pt_cnt[k]*COORD_DIM, &BB_let.pt_coord[0]+BB_let.pt_dsp[k]*COORD_DIM, false);
                                    spt_id[indx].ReInit(BB_let.pt_cnt[k]          , &BB_let.pt_id   [0]+BB_let.pt_dsp[k]          , false);
                                }
                            }
                            indx++;
                        }
            }

            { // Find near pairs
                std::vector<size_t> pair_vec;
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
                            if(r2<r2_near && spt_id[k][s]!=tpt_id[t]){
                                size_t pt_id_s=spt_id[k][s];
                                pair_vec.push_back(pt_id_s);
                            }
                        }
                    std::sort(&pair_vec[0],&pair_vec[0]+pair_vec.size());

                    for(size_t s=0;s<pair_vec.size();s++){
                        std::pair<size_t, size_t> new_pair;
                        new_pair.first=tdsp+t;
                        new_pair.second=pair_vec[s];
                        near_pair.push_back(new_pair);
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
    pvfmm::Vector<size_t> near_cnt(omp_p_);
    pvfmm::Vector<size_t> near_dsp(omp_p_); near_dsp[0]=0;
    for(size_t tid=0;tid<omp_p_;tid++){
        if(tid)
            near_dsp[tid]=near_pair_omp[tid-1].size()+near_dsp[tid-1];
        near_cnt[tid]=near_pair_omp[tid  ].size();
        near_size   +=near_pair_omp[tid  ].size();
    }

    near_trg_pt_id.ReInit(near_size);
    near_src_pt_id.ReInit(near_size);

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p_;tid++){
        size_t dsp=near_dsp[tid];
        size_t cnt=near_cnt[tid];
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];
        for(size_t i=0;i<cnt;i++){
            size_t loc_trg_id=near_pair[i].first;
            near_trg_pt_id[(dsp+i)]=pt_id[loc_trg_id];
            near_src_pt_id[(dsp+i)]=near_pair[i].second;
        }
    }
    pvfmm::Profile::Toc();
}
*/

template<typename Real_t>
void VesBoundingBox<Real_t>::FindNearPair(TREEGRID &BB_let, 
        pvfmm::Vector<size_t> &near_trg_box_id, pvfmm::Vector<size_t> &near_src_box_id)
{
    pvfmm::Profile::Tic("NearPair",&comm,true);
    PVFMMVec_t &box_min = BB_let.box_min;
    PVFMMVec_t &box_max = BB_let.box_max;
    pvfmm::Vector<size_t> &pt_id = BB_let.pt_id;
    pvfmm::Vector<size_t> &box_id = BB_let.box_id;
    pvfmm::Vector<pvfmm::MortonId> &tree_mid = BB_let.mid;
    pvfmm::Vector<size_t> &tree_pt_cnt = BB_let.pt_cnt;
    pvfmm::Vector<size_t> &tree_pt_dsp = BB_let.pt_dsp;

    std::vector<std::vector<std::pair<size_t, size_t> > > near_pair_omp(omp_p_); // (box_id, box_id)
    #pragma omp parallel num_threads(omp_p_)
    {
        size_t tid=omp_get_thread_num();
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];
        size_t tree_depth; Real_t r2_near;
        { // Set tree_depth, r_near
            tree_depth=BB_let.mid[0].GetDepth();
            r2_near=r_near_;
            r2_near*=r2_near;
        }
        Real_t s=pow(0.5,tree_depth);

        size_t FLOP=0;
        size_t a=((tid+0)*tree_mid.Dim())/omp_p_;
        size_t b=((tid+1)*tree_mid.Dim())/omp_p_;
        for(size_t i=a;i<b;i++){
            size_t tcnt=tree_pt_cnt[i];
            size_t tdsp=tree_pt_dsp[i];
            PVFMMVec_t tmin, tmax;
            pvfmm::Vector<size_t> tpt_id;
            pvfmm::Vector<size_t> tbox_id;
            if(tcnt){ // Set t_coord
                tmin.ReInit(tcnt*COORD_DIM,  &box_min[tdsp*COORD_DIM],false);
                tmax.ReInit(tcnt*COORD_DIM,  &box_max[tdsp*COORD_DIM],false);
                tpt_id.ReInit(tcnt,          &pt_id[tdsp]            ,false);
                tbox_id.ReInit(tcnt,         &box_id[tdsp]           ,false);
            }

            { // Find near pairs
                for(size_t t1=0;t1<tcnt;t1++){
                    for(size_t t2=0;t2<tcnt;t2++){
                        if( (tbox_id[t1]!=tbox_id[t2]) &&
                            CheckBBCollision(&tmin[t1*COORD_DIM], &tmax[t1*COORD_DIM], &tmin[t2*COORD_DIM],&tmax[t2*COORD_DIM])
                          )
                        {
                            std::pair<size_t, size_t> new_pair;
                            new_pair.first =  tbox_id[t1];
                            new_pair.second = tbox_id[t2];
                            near_pair.push_back(new_pair);
                        }
                    }

                    std::sort(near_pair.begin(),near_pair.end());
                    near_pair.erase(std::unique(near_pair.begin(), near_pair.end()), near_pair.end());
                
                    // TODO: add FLOP
                }
            }
        }
        pvfmm::Profile::Add_FLOP(FLOP);
    }

    size_t near_size=0;
    pvfmm::Vector<size_t> near_cnt(omp_p_);
    pvfmm::Vector<size_t> near_dsp(omp_p_); near_dsp[0]=0;
    for(size_t tid=0;tid<omp_p_;tid++){
        if(tid)
            near_dsp[tid]=near_pair_omp[tid-1].size()+near_dsp[tid-1];
        near_cnt[tid]=near_pair_omp[tid  ].size();
        near_size   +=near_pair_omp[tid  ].size();
    }

    near_trg_box_id.ReInit(near_size);
    near_src_box_id.ReInit(near_size);

    #pragma omp parallel for
    for(size_t tid=0;tid<omp_p_;tid++){
        size_t dsp=near_dsp[tid];
        size_t cnt=near_cnt[tid];
        std::vector<std::pair<size_t, size_t> >& near_pair=near_pair_omp[tid];
        for(size_t i=0;i<cnt;i++){
            near_trg_box_id[(dsp+i)]=near_pair[i].first;
            near_src_box_id[(dsp+i)]=near_pair[i].second;
        }
    }

    // scatter to box_ids
    size_t box_id_offset;
    {
        long long disp = 0;
        long long size = N_bbox_;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        box_id_offset = disp - size;
    }
    pvfmm::Vector<size_t> scatter_idx;
    pvfmm::par::SortScatterIndex(near_src_box_id, scatter_idx, comm, &box_id_offset);
    pvfmm::par::ScatterForward(near_src_box_id, scatter_idx, comm);
    pvfmm::par::ScatterForward(near_trg_box_id, scatter_idx, comm);

    pvfmm::Profile::Toc();
}

template<typename Real_t>
void VesBoundingBox<Real_t>::GlobalBoundingBox(Real_t *scale_xr, Real_t *shift_xr)
{
    Real_t& scale_x = *scale_xr;
    Real_t* shift_x = shift_xr;

    double loc_min_x[COORD_DIM];
    double loc_max_x[COORD_DIM];
    
    assert(N_bbox_>0);
    size_t n_src = N_bbox_;

    for(size_t k=0;k<COORD_DIM;k++){
        loc_min_x[k] = BB_min_[k];
        loc_max_x[k] = BB_max_[k];
    }
    
    for(size_t i=0;i<n_src;i++){
        const Real_t* x_min_=&BB_min_[i*COORD_DIM];
        const Real_t* x_max_=&BB_max_[i*COORD_DIM];
        for(size_t k=0;k<COORD_DIM;k++){
            if(loc_min_x[k]>x_min_[0]) loc_min_x[k]=x_min_[0];
            if(loc_max_x[k]<x_max_[0]) loc_max_x[k]=x_max_[0];
            ++x_min_; ++x_max_;
        }
    }

    double min_x[COORD_DIM];
    double max_x[COORD_DIM];
    MPI_Allreduce(loc_min_x, min_x, COORD_DIM, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(loc_max_x, max_x, COORD_DIM, MPI_DOUBLE, MPI_MAX, comm);

    static Real_t eps=machine_eps<Real_t>()*64; // Points should be well within the box.
    scale_x=1.0/(max_x[0]-min_x[0]+2*eps);
    for(size_t k=0;k<COORD_DIM;k++){
        scale_x=std::min(scale_x,(Real_t)(1.0/(max_x[k]-min_x[k]+2*eps)));
    }
    if(scale_x*0.0!=0.0) scale_x=1.0; // fix for scal_x=inf
    for(size_t k=0;k<COORD_DIM;k++){
        shift_x[k]=-min_x[k]*scale_x+eps;
    }
}
    
template<typename Real_t>
void VesBoundingBox<Real_t>::GenerateBBPoints()
{
    std::vector<size_t> bbox_nxyz(COORD_DIM*N_bbox_);
    std::vector<size_t> bbox_n(N_bbox_);
    std::vector<size_t> bbox_ndsp(N_bbox_);
    size_t n_sum = 0;

    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*N_bbox_)/omp_p_;
        size_t b = ((tid+1)*N_bbox_)/omp_p_;

        for(size_t i=a; i<b; i++){
            const Real_t* min_i = &BB_min_[i*COORD_DIM];
            const Real_t* max_i = &BB_max_[i*COORD_DIM];
            size_t* bbox_nxyzi     = &bbox_nxyz[i*COORD_DIM];

            int nx = std::ceil((max_i[0] - min_i[0])/leaf_size_);
            int ny = std::ceil((max_i[1] - min_i[1])/leaf_size_);
            int nz = std::ceil((max_i[2] - min_i[2])/leaf_size_);
            
            assert(nx > 0);assert(ny > 0);assert(nz > 0);
            
            bbox_nxyzi[0] = nx;
            bbox_nxyzi[1] = ny;
            bbox_nxyzi[2] = nz;
            bbox_n[i] = (nx+1)*(ny+1)*(nz+1);
        }
    }
        
    bbox_ndsp[0]=0; pvfmm::omp_par::scan(&bbox_n[0], &bbox_ndsp[0], bbox_n.size());
    n_sum = pvfmm::omp_par::reduce(&bbox_n[0], bbox_n.size());

    BB_pts_.ReInit(COORD_DIM*n_sum);
    BB_pts_min_.ReInit(COORD_DIM*n_sum);
    BB_pts_max_.ReInit(COORD_DIM*n_sum);
    BB_id_.ReInit(n_sum);

    size_t box_id_offset;
    {
        long long disp = 0;
        long long size = N_bbox_;
        MPI_Scan(&size, &disp, 1, MPI_LONG_LONG, MPI_SUM, comm);
        box_id_offset = disp - size;
    }

    #pragma omp parallel for
    for(size_t tid=0; tid<omp_p_; tid++){
        size_t a = ((tid+0)*N_bbox_)/omp_p_;
        size_t b = ((tid+1)*N_bbox_)/omp_p_;

        for(size_t i=a; i<b; i++){
            const Real_t* min_i = &BB_min_[i*COORD_DIM];
            const Real_t* max_i = &BB_max_[i*COORD_DIM];

            size_t box_id = box_id_offset + i;
            size_t nx = bbox_nxyz[i*COORD_DIM + 0];
            size_t ny = bbox_nxyz[i*COORD_DIM + 1];
            size_t nz = bbox_nxyz[i*COORD_DIM + 2];

            size_t ni = 0;
            size_t disp = bbox_ndsp[i];
            for(size_t nxi=0; nxi<nx+1; nxi++)
            for(size_t nyi=0; nyi<ny+1; nyi++)
            for(size_t nzi=0; nzi<nz+1; nzi++)
            {
                Real_t* BB_pts = &BB_pts_[COORD_DIM*(disp+ni)];
                BB_pts[0] = nxi*((max_i[0]-min_i[0])/nx) + min_i[0];
                BB_pts[1] = nyi*((max_i[1]-min_i[1])/ny) + min_i[1];
                BB_pts[2] = nzi*((max_i[2]-min_i[2])/nz) + min_i[2];

                Real_t* BB_pts_min = &BB_pts_min_[COORD_DIM*(disp+ni)];
                BB_pts_min[0] = min_i[0]; BB_pts_min[1] = min_i[1]; BB_pts_min[2] = min_i[2];
                Real_t* BB_pts_max = &BB_pts_max_[COORD_DIM*(disp+ni)];
                BB_pts_max[0] = max_i[0]; BB_pts_max[1] = max_i[1]; BB_pts_max[2] = max_i[2];
                
                BB_id_[disp+ni] = box_id;
                ni++;
            }
        }
    }
}

template<typename Real_t>
bool VesBoundingBox<Real_t>::CheckBBCollision(const Real_t *minA, const Real_t *maxA, const Real_t *minB, const Real_t *maxB)
{
    return ( (minA[0]<=maxB[0] && maxA[0]>= minB[0]) &&
             (minA[1]<=maxB[1] && maxA[1]>= minB[1]) &&
             (minA[2]<=maxB[2] && maxA[2]>= minB[2])
           );
}
