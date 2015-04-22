/**
 * @file   Scalars.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief The header file for the Scalars class. The implementation
 * of this class is in src/Scalars.cc.
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

#ifndef _SCALARS_H_
#define _SCALARS_H_

#include <utility>  // for pair
#include <cstddef>  // for size_t
#include "Array.h"
#include "Enums.h"
#include "Spharm.h"

/**
 * @class Scalar
 * @param T       Content datatype (float,double,etc.)
 * @param DT      Device type
 * @param DEVICE  Reference to the device object
 * @brief Container for scalar fields over spherical harmonics basis.
 *
 * @warning Defined iterators are not be use for direct iteration
 * because the data is on the device. They are intended for functions
 * that are aware of the device.
 */
template <typename T,    // content type
          typename DT,   // device class
          const DT &DEVICE>
class Scalars : public Array<T, DT, DEVICE>
{
  public:
    typedef Array<T, DT, DEVICE> array_type;
    typedef typename array_type::value_type     value_type;
    typedef typename array_type::device_type    device_type;
    typedef typename array_type::iterator       iterator;
    typedef typename array_type::const_iterator const_iterator;

    //! if grid_dim is EMPTY_GRID, it is defaulted to the implied
    //! value throgh sh_
    explicit Scalars(size_t num_subs = 0, int sh_order = -1,
        std::pair<int, int> grid_dim = EMPTY_GRID);

    explicit Scalars(std::istream &is, Streamable::Format format);

    virtual ~Scalars();

    //! The number of dimensions of the field
    static inline int getTheDim();

    //! Spherical harmonics order
    inline int getShOrder() const;

    //! The u-v grid dimension size
    inline std::pair<int, int> getGridDim() const;

    //! Stride between sub-fields
    inline size_t getStride() const;

    //! number of sub-components
    virtual inline size_t getNumSubs() const;

    //! number of sub-fields each a scalar function
    inline size_t getNumSubFuncs() const;

    //! Total length of data for each sub-component
    virtual inline size_t getSubLength() const;

    //! Total length of data for each sub-field
    inline size_t getSubFuncLength() const;

    //! required size of the underlying array
    inline size_t req_arr_size() const;

    //! Resizing the container. This doesn't interpolate
    inline void resize(size_t new_num_subs,
        int new_sh_order = -1,
        std::pair<int, int> new_grid_dim = EMPTY_GRID);

    //! resizing so that this is compatible with rhs
    inline void replicate(Scalars<T, DT, DEVICE> const& rhs);

    //! resizing so that this is the same size as rhs
    inline void match_size(Scalars<T, DT, DEVICE> const& rhs);

    //! head of each sub component
    inline iterator getSubN_begin(size_t n);
    //! tail of each sub component
    inline iterator getSubN_end(size_t n);

    //! head of each sub function
    inline iterator getSubFuncN_begin(size_t n);
    //! head of each sub function
    inline iterator getSubFuncN_end(size_t n);

    //! head of each sub component
    inline const_iterator getSubN_begin(size_t n) const;
    //! tail of each sub component
    inline const_iterator getSubN_end(size_t n) const;

    //! head of each sub fuction
    inline const_iterator getSubFuncN_begin(size_t n) const;
    //! tail of each sub function
    inline const_iterator getSubFuncN_end(size_t n) const;

    // From streamable class --------------------------------------------------
    // ------------------------------------------------------------------------
    virtual Error_t pack(std::ostream &os, Streamable::Format format) const;
    virtual Error_t unpack(std::istream &is, Streamable::Format format);

  protected:
    //! set number of sub-components
    virtual inline void setNumSubs(size_t num_subs);

    int sh_order_;
    std::pair<int, int> grid_dim_;
    size_t stride_;
    size_t num_sub_funcs_;

    static const int scalar_dim_ = 1; //scalar field has dim 1

    //! protected to copy constructor to enforce pass by ref
    Scalars(Scalars<T, DT, DEVICE> const& sc_in);

    //! protected to assignment to enforce pass by ref
    Scalars<T, DT, DEVICE>& operator=(const Scalars<T, DT, DEVICE>& rhs);
};

//! location of flag for streaming
static int scalars_xalloc(std::ios_base::xalloc());

//! Streaming manipulator
std::ios_base& toggle_show_entries(std::ios_base& os)
{
    os.iword(scalars_xalloc) = !os.iword(scalars_xalloc);
    return os;
}

//! Utitly for streaming
template <typename T, typename DT, const DT &DEVICE>
std::ostream& operator<<(std::ostream& output,
    const Scalars<T, DT, DEVICE>&sc);

#include "Scalars.cc"

#endif //_SCALARS_H_
