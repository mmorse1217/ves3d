/**
 * @file   Scalars.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 13:48:19 2010
 * 
 * @brief  Implementation for the Scalars class. 
 */

template<typename T, enum DeviceType DT> 
Scalars<T, DT>::Scalars() :
    device_(NULL),
    grid_dim_(EMPTY_GRID)
    data_(NULL)
{}

template<typename T, enum DeviceType DT> 
Scalars<T, DT>::Scalars(Device<DT> *device, int sh_order, 
    size_t num_funs, pair<int, int> grid_dim) :
    device_(device),
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? GridDimOf(sh_order_) : grid_dim),
    fun_length_(grid_dim_.first * grid_dim_.second),
    num_funs_(num_funs),
    capacity_(fun_length_ * num_funs_),
    data_((T*) device_->Malloc(the_dim_ * capacity_ * sizeof(T)))
{}

template<typename T, enum DeviceType DT> 
Scalars<T, DT>::~Scalars()
{
    device_->Free(data_);
}

template<typename T, enum DeviceType DT> 
void Scalars<T,DT>::Resize(size_t new_num_funs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    sh_order_ = (new_sh_order == 0) ? sh_order_ : new_sh_order;
    grid_dim_ = (new_grid_dim == EMPTY_GRID) ? GridDimOf(sh_order_) : new_grid_dim;
    
    size_t new_stride(grid_dim_.first * grid_dim_.second);
    size_t new_capacity(new_stride * new_num_funs);
    
    if(new_capacity > capacity_)
    {
        T *data_old(data_);
        data_ = (T*) device_->Malloc(the_dim_ * new_capacity * sizeof(T));

        if(data_old != NULL) 
            device_->Memcpy(data_, data_old, 
                the_dim_ * fun_length_ * num_funs_ * sizeof(T), 
                MemcpyDeviceToDevice);
        device_->Free(data_old);
        capacity_ = new_capacity;
    }    

    fun_length_ = new_stride;
    num_funs_ = new_num_funs;
}

template<typename T, enum DeviceType DT> 
const Device<DT>* Scalars<T,DT>::GetDevicePtr() const
{
    return(device_);
}

template<typename T, enum DeviceType DT> 
int Scalars<T,DT>::GetShOrder() const
{
    return(sh_order_);
}

template<typename T, enum DeviceType DT> 
pair<int, int> Scalars<T,DT>::GetGridDim() const
{
    return(grid_dim_);
}

template<typename T, enum DeviceType DT> 
size_t Scalars<T,DT>::GetFunLength() const
{
    return(fun_length_);
}

template<typename T, enum DeviceType DT> 
size_t Scalars<T,DT>::GetNumFuns() const
{
    return(num_funs_);
}

template<typename T, enum DeviceType DT> 
size_t Scalars<T,DT>::Size() const
{
    return(the_dim_ * fun_length_ * num_funs_);
}

template<typename T, enum DeviceType DT> 
T& Scalars<T,DT>::operator[](size_t idx)
{
    return(data_[idx]);
}

template<typename T, enum DeviceType DT> 
const T& Scalars<T,DT>::operator[](size_t idx) const
{
    return(data_[idx]);
}


template<typename T, enum DeviceType DT> 
Scalars<T,DT>::iterator Scalars<T,DT>::begin()
{
    return(data_);
}

template<typename T, enum DeviceType DT> 
Scalars<T,DT>::const_iterator Scalars<T,DT>::begin() const
{
    return(data_);
}

template<typename T, enum DeviceType DT> 
Scalars<T,DT>::iterator Scalars<T,DT>::end()
{
    return(data_ + the_dim_ * fun_length_ * num_funs_);
}

template<typename T, enum DeviceType DT> 
Scalars<T,DT>::const_iterator Scalars<T,DT>::end() const
{
    return(data_ + the_dim_ * fun_length_ * num_funs_);
}
