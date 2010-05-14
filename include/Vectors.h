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

template <typename T, enum DeviceType DT> 
class Vectors : public Scalars<T,DT>
{
  protected:
    size_t num_vecs_;
    enum CoordinateOrder point_order_;
    static const int the_dim_ = 3;
    
  public: 
    Vectors();
    Vectors(Device<DT> *device, int sh_order = 0,
        size_t num_vecs = 0, pair<int, int> grid_dim = EMPTY_GRID);
    
    virtual void Resize(size_t new_num_vecs, int new_sh_order = 0,
        pair<int, int> new_grid_dim = EMPTY_GRID);

    inline size_t GetNumVecs() const;
    inline void SetPointOrder(enum CoordinateOrder new_order);    
    inline enum CoordinateOrder GetPointOrder() const;
};

#include "Vectors.cc"
#endif
