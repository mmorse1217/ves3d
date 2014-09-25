/**
 * @file   Vector.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief Container of vector fields of funtions over spherical
 * harmonic grid
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _VECTORS_H_
#define _VECTORS_H_

#include "Scalars.h"
#include "Enums.h"

template <typename T,    // content type
          typename DT,   // device class
          const DT &DEVICE>
class Vectors : public Scalars<T, DT, DEVICE>
{
  public:
    typedef Scalars<T, DT, DEVICE> scalars_type;
    typedef typename scalars_type::value_type     value_type;
    typedef typename scalars_type::device_type    device_type;
    typedef typename scalars_type::iterator       iterator;
    typedef typename scalars_type::const_iterator const_iterator;

    explicit Vectors(
        size_t num_vecs = 0,
        int sh_order = -1,
        std::pair<int, int> grid_dim = EMPTY_GRID,
        CoordinateOrder po = AxisMajor);

    //! The number of dimensions of the field
    static inline int getTheDim();

    //! number of sub vector fields
    virtual inline size_t getNumSubs() const;

    //! Total length of data for each sub-field
    virtual inline size_t getSubLength() const;

    //! The ordering of vectors when iterated upon
    //! Point ordering could be AxisMajor or PointMajor
    //! @todo call to this should either transpose the points or store
    //! the info and iterate on points with stride
    inline void setPointOrder(enum CoordinateOrder new_order);
    inline enum CoordinateOrder getPointOrder() const;

  protected:
    //! set number of sub-components
    virtual inline void setNumSubs(size_t num_subs);

    enum CoordinateOrder point_order_;
    static const int vector_dim_ = 3;
};

#include "Vectors.cc"

#endif
