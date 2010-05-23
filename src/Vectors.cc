/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */

template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>  
Vectors<T, DT, DEVICE>::Vectors(int sh_order,
    size_t num_vecs, pair<int, int> grid_dim) :
    Scalars<T, DT, DEVICE>(sh_order, num_vecs * the_dim_, grid_dim), 
    num_vecs_(num_vecs),
    point_order_(AxisMajor)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::Resize(size_t new_num_vecs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    Scalars<T, DT, DEVICE>::Resize(the_dim_ * new_num_vecs, 
        new_sh_order, new_grid_dim);
    num_vecs_ = new_num_vecs;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Vectors<T, DT, DEVICE>::GetNumVecs() const
{
    return(num_vecs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::SetPointOrder(enum CoordinateOrder new_order)
{
    point_order_ = new_order;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>  
enum CoordinateOrder Vectors<T, DT, DEVICE>::GetPointOrder() const
{
    return(point_order_);
}


