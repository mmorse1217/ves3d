/**
 * @file   Vectors.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Fri Jan 29 09:36:58 2010
 * 
 * @brief  The implementation of Vectors.
 */

template <typename T, enum DeviceType DT, const Device<DT> &DEVICE>  
Vectors<T, DT, DEVICE>::Vectors(size_t num_vecs, int sh_order,
    pair<int, int> grid_dim) :
    Scalars<T, DT, DEVICE>(num_vecs * this->the_dim_, sh_order, grid_dim), 
    num_vecs_(num_vecs),
    point_order_(AxisMajor)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::resize(size_t new_num_vecs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    Scalars<T, DT, DEVICE>::resize(this->the_dim_ * new_num_vecs, 
        new_sh_order, new_grid_dim);
    num_vecs_ = new_num_vecs;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Vectors<T, DT, DEVICE>::getNumSubs() const
{
    return(this->num_vecs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Vectors<T, DT, DEVICE>::getTheDim() const
{
    return(this->the_dim_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Vectors<T, DT, DEVICE>::replicate(Vectors<T, DT, DEVICE> const& vec_in)
{
    this->resize(vec_in.getNumSubs(), vec_in.getShOrder(), 
        vec_in.getGridDim());
    
    point_order_ = vec_in.getPointOrder();    
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


