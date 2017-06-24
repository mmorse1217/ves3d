#ifndef _VES_BOUNDING_BOX_H_
#define _VES_BOUNDING_BOX_H_
//#define VES_BOUNDING_BOX_DEBUG

#include <mpi.h>
#include <vector.hpp>
#include <matrix.hpp>

template<typename Real_t>
class VesBoundingBox{
    public:
        typedef typename pvfmm::Vector<Real_t> PVFMMVec_t;
        // a tree <=> a grid
        typedef struct TREEGRID {
            pvfmm::Vector<pvfmm::MortonId> mid;   // mortonid of leaf nodes
            pvfmm::Vector<size_t> pt_cnt;         // point count
            pvfmm::Vector<size_t> pt_dsp;         // point displ
            pvfmm::Vector<pvfmm::MortonId> mins;  // first non-ghost node

            pvfmm::Vector<size_t> pt_id;          // scatter id of points
            pvfmm::Vector<size_t> box_id;         // box id for each point
            PVFMMVec_t box_min;                   // mins
            PVFMMVec_t box_max;                   // maxs
        } TREEGRID;
 
        // constructor, set box_size_, comm
        VesBoundingBox(Real_t box_size=-1, MPI_Comm c=MPI_COMM_WORLD);

        // set bounding boxes for vesicles, the bounding boxes are space-time bounding box with minimum separation distance
        void SetVesBoundingBox(const PVFMMVec_t& ves_coord_s, const PVFMMVec_t& ves_coord_e, 
                const Real_t min_sep, const int sh_order, const int sh_order_up);
        template<typename Vec>
        void SetVesBoundingBox(const Vec& ves_coord_s, const Vec& ves_coord_e, 
                const Real_t min_sep, const int sh_order, const int sh_order_up);
        // set bounding boxes provided with two points for each bounding box: BB_min, BB_max
        void SetVesBoundingBox(const PVFMMVec_t& BB_min, const PVFMMVec_t& BB_max);
        
        // get all pairs of intersecting bounding boxes
        void GetContactBoundingBoxPair(std::vector< std::pair<size_t, size_t> > &BBIPairs);
        
        // set the grid parameters, bbox_(shift,scale), r_near_(intrinsic size of grid box),
        // leaf_size_  <=> tree_depth_ (grid size) TODO: eliminate either leaf_size_ or tree_depth_
        void SetTreeParams(); // set bbox_, r_near_, tree_depth_, leaf_size_
        
        // sort the generated points to processes, construct a local grid, each grid box is represented by a morton id
        void ConstructLocalTree(TREEGRID &BB_let);
        
        // bring the required ghost grid boxex
        void AddGhostNodes(TREEGRID &BB_let);

        // use the constructed grid to find pairs(box_id_i, box_id_j) of intersecting bounding boxes
        void FindNearPair(TREEGRID &BB_let, std::vector< std::pair<size_t, size_t> > &BBIPairs);

        // find the global bounding box of all the bounding boxes; find the shift and scale so that each point can be
        // transferred to a unit box
        void GlobalBoundingBox(Real_t *scale_xr, Real_t *shift_xr);
        
        // generate points for each bounding box, the distance between two adjacent points of a bounding box 
        // should be less than or equal the leaf_size(grid box size)
        void GenerateBBPoints();

        // check if two bounding boxes intersect
        bool CheckBBCollision(const Real_t *minA, const Real_t *maxA, const Real_t *minB, const Real_t *maxB);

    private:
        // TODO: some constructor, assignment operator
        //VesBoundingBox(const VesBoundingBox &);
        //VesBoundingBox& operator=(const VesBoundingBox &);

        Real_t box_size_; // store the periodic box_size
        MPI_Comm comm; // communicator of MPI

        // the following three containers are of size COORD_DIM*N_pts
        // N_pts is the total number of points generated of all bounding boxes
        PVFMMVec_t BB_pts_; // all generated points of all bounding boxes
        // there are duplicates for points of same bounding box, TODO: is there a compact representation without duplicates
        PVFMMVec_t BB_pts_min_; // for each generated point, store the min of its belonging bounding box
        PVFMMVec_t BB_pts_max_; // for each generated point, store the max of its belonging bounding box

        // the following container are of size N_pts
        // there are duplicates for points of same bounding box, TODO: is there a compact representation without duplicates
        pvfmm::Vector<size_t> BB_id_; // for each generated point, store the bounding box id of its belonging bounding box
        
        // bounding boxes, each bounding box has a min_pt and a max_pt
        // the following two containers are of size N_bbox_
        PVFMMVec_t BB_min_;
        PVFMMVec_t BB_max_;
    
        Real_t bbox_[4]; // {x,y,z,s} : shift, scale
        // TODO: choose a "good" size so that the assumption(each process should have at least one (morton id) grid box) 
        // won't fail, or remove the assumption
        Real_t r_near_; // the intrinsic size of grid box, r_near_ = bb_size_max or bb_size_min or bb_size_avg
        Real_t leaf_size_; // grid box size, for free space, leaf_size_ = r_near_/2; for periodic case, r_near_/4 < leaf_size_ <= r_near/2
        size_t tree_depth_; // depth of the "tree"/grid, level 0 - a grid of one grid box, level n - a grid of 2^n grid boxes

        // number of MPI processes and current process id
        int np_, rank_;
        // number of total openmp threads
        size_t omp_p_;
        // number of bounding boxes
        size_t N_bbox_;
        // number of points
        size_t N_pts_;
};

#include "VesBoundingBox.cc"

#endif // _VES_BOUNDING_BOX_H_
