/**
 * @file   Scalars.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 13:48:19 2010
 * 
 * @brief  Implementation for the Scalars class. 
 */

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::Scalars(int sh_order, size_t num_funs, 
    pair<int, int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? GridDimOf(sh_order_) : grid_dim),
    fun_length_(grid_dim_.first * grid_dim_.second),
    num_funs_(num_funs),
    capacity_(fun_length_ * num_funs_),
    data_((T*) DEVICE.Malloc(the_dim_ * capacity_ * sizeof(T)))
{}

template<typename T, enum DeviceType DT,
const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::Scalars(Scalars<T, DT, DEVICE> const& sc_in) :
    sh_order_(sc_in.sh_order_),
    grid_dim_(sc_in.grid_dim_),
    fun_length_(grid_dim_.first * grid_dim_.second),
    num_funs_(sc_in.num_funs_),
    capacity_(fun_length_ * num_funs_),
    data_((T*) DEVICE.Malloc(the_dim_ * capacity_ * sizeof(T)))
{
    DEVICE.Memcpy(data_, sc_in.begin(), sc_in.Size() * sizeof(T), 
        MemcpyDeviceToDevice);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::~Scalars()
{
    DEVICE.Free(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Scalars<T, DT, DEVICE>::Resize(size_t new_num_funs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    sh_order_ = (new_sh_order == 0) ? sh_order_ : new_sh_order;
    grid_dim_ = (new_grid_dim == EMPTY_GRID) ? GridDimOf(sh_order_) : new_grid_dim;
    
    size_t new_stride(grid_dim_.first * grid_dim_.second);
    size_t new_capacity(new_stride * new_num_funs);
    
    if(new_capacity > capacity_)
    {
        T *data_old(data_);
        data_ = (T*) DEVICE.Malloc(the_dim_ * new_capacity * sizeof(T));

        if(data_old != NULL) 
            DEVICE.Memcpy(data_, data_old, 
                the_dim_ * fun_length_ * num_funs_ * sizeof(T), 
                MemcpyDeviceToDevice);
        DEVICE.Free(data_old);
        capacity_ = new_capacity;
    }    

    fun_length_ = new_stride;
    num_funs_ = new_num_funs;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
const Device<DT>& Scalars<T, DT, DEVICE>::GetDevice() const
{
    return(DEVICE);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Scalars<T, DT, DEVICE>::GetShOrder() const
{
    return(sh_order_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
pair<int, int> Scalars<T, DT, DEVICE>::GetGridDim() const
{
    return(grid_dim_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::GetFunLength() const
{
    return(fun_length_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::GetNumFuns() const
{
    return(num_funs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::Size() const
{
    return(the_dim_ * fun_length_ * num_funs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
T& Scalars<T, DT, DEVICE>::operator[](size_t idx)
{
    return(data_[idx]);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
const T& Scalars<T, DT, DEVICE>::operator[](size_t idx) const
{
    return(data_[idx]);
}


template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::iterator Scalars<T, DT, DEVICE>::begin()
{
    return(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::const_iterator Scalars<T, DT, DEVICE>::begin() const
{
    return(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::iterator Scalars<T, DT, DEVICE>::end()
{
    return(data_ + the_dim_ * fun_length_ * num_funs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::const_iterator Scalars<T, DT, DEVICE>::end() const
{
    return(data_ + the_dim_ * fun_length_ * num_funs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
std::ostream& operator<<(std::ostream& output, Scalars<T, DT, DEVICE> &sc)
{
    output<<"=====================================================\n"
          <<"SH order            : "<<sc.GetShOrder()<<"\n"
          <<"Grid size           : ( "
          <<sc.GetGridDim().first<<" , "<<sc.GetGridDim().second<<" )"<<"\n"
          <<"stride              : "<<sc.GetFunLength()<<"\n"
          <<"Number of functions : "<<sc.GetNumFuns()<<"\n"
          <<"=====================================================\n";
    for(typename Scalars<T,DT,DEVICE>::iterator it = sc.begin(); it !=sc.end(); ++it)
    {
        output<<*it<<" ";

        if((it-sc.begin() + 1)%sc.GetGridDim().second == 0)
            output<<endl;
        
        if((it-sc.begin() + 1)%sc.GetFunLength() == 0)
            output<<endl;

    }
    return(output);
}
