/**
 * @file   Vector.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
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

template <typename T, typename DT, const DT &DEVICE>
Vectors<T, DT, DEVICE>::Vectors(size_t num_sub_vecs, int sh_order,
    std::pair<int, int> grid_dim, CoordinateOrder po) :
    Scalars<T, DT, DEVICE>(num_sub_vecs * this->the_dim_,
                           sh_order,
                           grid_dim),
    point_order_(po)
{}

template <typename T, typename DT, const DT &DEVICE>
int Vectors<T, DT, DEVICE>::getTheDim()
{
    return(the_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Vectors<T, DT, DEVICE>::getNumSubs() const
{
    return(scalars_type::num_subs_ / the_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Vectors<T, DT, DEVICE>::getSubLength() const
{
    return((scalars_type::num_subs_ > 0 ) * the_dim_ * this->getStride() );
}

template <typename T, typename DT, const DT &DEVICE>
void Vectors<T, DT, DEVICE>::resize(size_t new_num_sub_vecs,
    int new_sh_order,
    std::pair<int, int> new_grid_dim)
{
    Scalars<T, DT, DEVICE>::resize(this->getTheDim() * new_num_sub_vecs,
        new_sh_order, new_grid_dim);
}

template <typename T, typename DT, const DT &DEVICE>
void Vectors<T, DT, DEVICE>::replicate(
    Vectors<T, DT, DEVICE> const& other)
{
    this->resize(other.getNumSubs(), other.getShOrder(),
        other.getGridDim());
    this->setPointOrder(other.getPointOrder());
}

template <typename T, typename DT, const DT &DEVICE>
void Vectors<T, DT, DEVICE>::setPointOrder(enum CoordinateOrder new_order)
{
    point_order_ = new_order;
}

template <typename T, typename DT, const DT &DEVICE>
enum CoordinateOrder Vectors<T, DT, DEVICE>::getPointOrder() const
{
    return(point_order_);
}

// template <typename T, typename DT, const DT &DEVICE>
// typename Vectors<T, DT, DEVICE>::iterator
// Vectors<T, DT, DEVICE>::getSubN_begin(size_t n)
// {
//     assert(n<num_sub_vecs_);
//     return( Scalars<T, DT, DEVICE>::getSubN_begin(size_t n*the_dim_));

//     return(this->begin() + n * this->getTheDim() * this->getStride());
// }

// template <typename T, typename DT, const DT &DEVICE>
// typename Vectors<T, DT, DEVICE>::const_iterator Vectors<T, DT, DEVICE>::getSubN(size_t n) const
// {
//     assert(n<num_sub_vecs_);
//     return(this->begin() + n * this->getTheDim() * this->getStride());
// }
