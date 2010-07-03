/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */
template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>  
Vectors<T, DT, DEVICE>::Vectors(size_t num_sub_vecs, int sh_order,
    pair<int, int> grid_dim) :
    Scalars<T, DT, DEVICE>(num_sub_vecs * this->the_dim_, sh_order, grid_dim), 
    num_sub_vecs_(num_sub_vecs),
    point_order_(AxisMajor)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::resize(size_t new_num_sub_vecs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    Scalars<T, DT, DEVICE>::resize(this->getTheDim() * new_num_sub_vecs, 
        new_sh_order, new_grid_dim);
    num_sub_vecs_ = new_num_sub_vecs;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Vectors<T, DT, DEVICE>::getNumSubs() const
{
    return(this->num_sub_vecs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Vectors<T, DT, DEVICE>::getTheDim()
{
    return(the_dim_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::replicate(Vectors<T, DT, DEVICE> const& vec_in)
{
    this->resize(vec_in.getNumSubs(), vec_in.getShOrder(), 
        vec_in.getGridDim());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::setPointOrder(enum CoordinateOrder new_order)
{
    ///@todo perform the transpose here.
    point_order_ = new_order;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>  
enum CoordinateOrder Vectors<T, DT, DEVICE>::getPointOrder() const
{
    return(point_order_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Vectors<T, DT, DEVICE>::iterator Vectors<T, DT, DEVICE>::getSubN(size_t n)
{
    assert(n<num_sub_vecs_);
    return(this->begin() + n * this->getTheDim() * this->getStride());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
typename Vectors<T, DT, DEVICE>::const_iterator Vectors<T, DT, DEVICE>::getSubN(size_t n) const
{
    assert(n<num_sub_vecs_);
    return(this->begin() + n * this->getTheDim() * this->getStride());
}
