/**
 * @file   Scalars.cc
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief  Implementation for the Scalars class.
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


inline std::pair<int, int> gridDimOf(int sh_order)
{
    return((sh_order >= 0) ?
        std::make_pair(sh_order + 1, 2 * sh_order + 2) :
        EMPTY_GRID);
}

template <typename T, typename DT, const DT &DEVICE>
Scalars<T, DT, DEVICE>::Scalars(size_t num_subs, int sh_order,
    std::pair<int, int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? gridDimOf(sh_order_) : grid_dim),
    stride_(grid_dim_.first * grid_dim_.second),
    num_subs_(num_subs)
{
    Array<T, DT, DEVICE>::resize(this->arr_size());
}

template <typename T, typename DT, const DT &DEVICE>
Scalars<T, DT, DEVICE>::~Scalars()
{}

template <typename T, typename DT, const DT &DEVICE>
int Scalars<T, DT, DEVICE>::getTheDim()
{
    return(the_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
int Scalars<T, DT, DEVICE>::getShOrder() const
{
    return(sh_order_);
}

template <typename T, typename DT, const DT &DEVICE>
std::pair<int, int> Scalars<T, DT, DEVICE>::getGridDim() const
{
    return(grid_dim_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Scalars<T, DT, DEVICE>::getStride() const
{
    return(stride_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Scalars<T, DT, DEVICE>::getNumSubs() const
{
    return(num_subs_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Scalars<T, DT, DEVICE>::getSubLength() const
{
    return((num_subs_ > 0 ) * the_dim_ * stride_);
}

template <typename T, typename DT, const DT &DEVICE>
size_t Scalars<T, DT, DEVICE>::arr_size() const
{
    return( the_dim_ * stride_ * num_subs_ );
}

template <typename T, typename DT, const DT &DEVICE>
void Scalars<T, DT, DEVICE>::resize(size_t new_num_subs,
    int new_sh_order,
    std::pair<int, int> new_grid_dim)
{
    num_subs_ = new_num_subs;
    sh_order_ = (new_sh_order == -1) ? sh_order_ : new_sh_order;
    grid_dim_ = (new_grid_dim == EMPTY_GRID) ?
        gridDimOf(sh_order_) : new_grid_dim;
    stride_   = grid_dim_.first * grid_dim_.second;

    Array<T, DT, DEVICE>::resize(this->arr_size());
}

template <typename T, typename DT, const DT &DEVICE>
void Scalars<T, DT, DEVICE>::replicate(Scalars<T, DT, DEVICE> const& sc_in)
{
    this->resize(sc_in.getNumSubs(), sc_in.getShOrder(),
        sc_in.getGridDim());
}

template <typename T, typename DT, const DT &DEVICE>
typename Scalars<T, DT, DEVICE>::iterator
Scalars<T, DT, DEVICE>::getSubN_begin(size_t n)
{
    return(Array<T, DT, DEVICE>::begin() +
        n * this->getTheDim() * this->getStride());
}

template <typename T, typename DT, const DT &DEVICE>
typename Scalars<T, DT, DEVICE>::iterator
Scalars<T, DT, DEVICE>::getSubN_end(size_t n)
{
    return(Array<T, DT, DEVICE>::begin() +
        (n+1) * this->getTheDim() * this->getStride());
}

template <typename T, typename DT, const DT &DEVICE>
typename Scalars<T, DT, DEVICE>::const_iterator
Scalars<T, DT, DEVICE>::getSubN_begin(size_t n) const
{
    return(Array<T, DT, DEVICE>::begin() +
        n * this->getTheDim() * this->getStride());
}

template <typename T, typename DT, const DT &DEVICE>
typename Scalars<T, DT, DEVICE>::const_iterator
Scalars<T, DT, DEVICE>::getSubN_end(size_t n) const
{
    return(Array<T, DT, DEVICE>::begin() +
        (n+1) * this->getTheDim() * this->getStride());
}

/////////////////////////////////////////////////////////////////////
template <typename T, typename DT, const DT &DEVICE>
std::ostream& operator<<(std::ostream& output, const Scalars<T, DT, DEVICE> &sc)
{
    output<<" =====================================================\n"
          <<"  SH order            : "<<sc.getShOrder()<<"\n"
          <<"  Grid size           : ( "
          <<sc.getGridDim().first<<" , "<<sc.getGridDim().second<<" )"<<"\n"
          <<"  stride              : "<<sc.getStride()<<"\n"
          <<"  Number of functions : "<<sc.getNumSubs()<<"\n"
          <<" =====================================================\n";

    return(output);
}

template <typename Container>
ShowEntries<Container>::ShowEntries(const Container &c) : c_(c) {}

template <typename Container>
std::ostream& ShowEntries<Container>::operator()(std::ostream &output) const
{
    output<<c_;

    typedef typename Container::value_type T;
    T *buffer(new T[this->c_.size()]);
    this->c_.getDevice().Memcpy(buffer,
        this->c_.begin(),
        this->c_.size() * sizeof(T),
        Container::device_type::MemcpyDeviceToHost);

    for(size_t ii=0; ii<this->c_.size(); ++ii)
    {
        output<<buffer[ii]<<" ";

        if((ii + 1)%this->c_.getGridDim().second == 0)
            output<<std::endl;

        if((ii + 1)%this->c_.getStride() == 0)
            output<<std::endl;
    }

    delete[] buffer;

    return(output);
}

template<typename Container>
std::ostream& operator<<(std::ostream& output, const ShowEntries<Container> &se)
{
    return(se(output));
}
