/**
 * @file   SHVectors.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Thu Jan 28 14:23:48 2010
 * 
 * @brief The header file for the Vectors class. See
 * src/Vectors.cc for its implementation.
 */

#ifndef _VECTORS_H_
#define _VECTORS_H_

#include "Scalars.h"

namespace containers
{
    template <typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
    class Vectors : public Scalars<T, DT, DEVICE>
    {
      protected:
        size_t num_vecs_;
        enum CoordinateOrder point_order_;
        
        static const int the_dim_ = 3;

      public: 
        Vectors(size_t num_vecs = 0, int sh_order = 0,
            pair<int, int> grid_dim = EMPTY_GRID);
        
        virtual void resize(size_t new_num_vecs, int new_sh_order = 0,
            pair<int, int> new_grid_dim = EMPTY_GRID);
        
        ///This is intentionally non-virtual
        inline size_t getNumSubs() const; 
        inline int getTheDim() const;

        inline void replicate(Vectors<T, DT, DEVICE> const& vec_in);        
        inline void setPointOrder(enum CoordinateOrder new_order);    
        inline enum CoordinateOrder getPointOrder() const;
    };

#include "Vectors.cc"
}

#endif
