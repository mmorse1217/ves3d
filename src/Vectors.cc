/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */

template <typename T, enum DeviceType DT> 
Vectors<T,DT>::Vectors() :
    num_vecs_(0),
    point_order_(AxisMajor)
{}

template <typename T, enum DeviceType DT> 
Vectors<T,DT>::Vectors(Device<DT> *device, int sh_order,
    size_t num_vecs, pair<int, int> grid_dim) :
    Scalars<T,DT>(device, sh_order, num_vecs * the_dim_, grid_dim),
    num_vecs_(num_vecs),
    point_order_(AxisMajor)
{}

template<typename T, enum DeviceType DT> 
void Vectors<T,DT>::Resize(size_t new_num_vecs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    Scalars<T,DT>::Resize(the_dim_ * new_num_vecs, new_sh_order, new_grid_dim);
    num_vecs_ = new_num_vecs;
}

template<typename T, enum DeviceType DT> 
size_t Vectors<T,DT>::GetNumVecs() const
{
    return(num_vecs_);
}

template<typename T, enum DeviceType DT> 
void Vectors<T,DT>::SetPointOrder(enum CoordinateOrder new_order)
{
    point_order_ = new_order;
}

template<typename T, enum DeviceType DT> 
enum CoordinateOrder Vectors<T,DT>::GetPointOrder() const
{
    return(point_order_);
}


