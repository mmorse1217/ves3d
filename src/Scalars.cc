/**
 * @file   Scalars.cc
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 26 13:48:19 2010
 * 
 * @brief  Implementation for the Scalars class. 
 */
template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::Scalars(size_t num_subs, int sh_order, 
    pair<int, int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? gridDimOf(sh_order_) : grid_dim),
    stride_(grid_dim_.first * grid_dim_.second),
    num_subs_(num_subs),
    capacity_(stride_ * num_subs_),
    data_((capacity_ > 0) ? (T*) DEVICE.Malloc(the_dim_ * capacity_ * sizeof(T)) : NULL)
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::~Scalars()
{
    DEVICE.Free(data_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Scalars<T, DT, DEVICE>::resize(size_t new_num_subs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    ///Note that num_subs = 0 frees all the memory (no recovery).
    
    sh_order_ = (new_sh_order == -1) ? sh_order_ : new_sh_order;
    grid_dim_ = (new_grid_dim == EMPTY_GRID) ? 
        gridDimOf(sh_order_) : new_grid_dim;
    
    size_t new_stride(grid_dim_.first * grid_dim_.second);
    size_t new_capacity(new_stride * new_num_subs);

    if ( new_capacity > capacity_ )
    {
        cout<<new_capacity<<" "<<capacity_<<endl;

        T *data_old(data_);
        data_ = (T*) DEVICE.Malloc(the_dim_ * new_capacity * sizeof(T));
        
        if(data_old != NULL) 
        {
            DEVICE.Memcpy(data_, data_old, 
                the_dim_ * stride_ * num_subs_ * sizeof(T), 
                MemcpyDeviceToDevice);
            DEVICE.Free(data_old);
        }
        capacity_ = new_capacity;
    } 
    else if (new_capacity == 0 && data_ != NULL)
    {
        DEVICE.Free(data_);
        data_ = NULL;
        capacity_ = 0;
    }

    stride_ = new_stride;
    num_subs_ = new_num_subs;
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
const Device<DT>& Scalars<T, DT, DEVICE>::getDevice()
{
    return(DEVICE);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Scalars<T, DT, DEVICE>::getShOrder() const
{
    return(sh_order_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
pair<int, int> Scalars<T, DT, DEVICE>::getGridDim() const
{
    return(grid_dim_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::getStride() const
{
    return(stride_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::size() const
{
    return(this->getTheDim() * this->getStride() * this->getNumSubs());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::getNumSubs() const
{
    return(num_subs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Scalars<T, DT, DEVICE>::getTheDim()
{
    return(the_dim_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Scalars<T, DT, DEVICE>::replicate(Scalars<T, DT, DEVICE> const& sc_in)
{
    this->resize(sc_in.getNumSubs(), sc_in.getShOrder(), 
        sc_in.getGridDim());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Scalars<T, DT, DEVICE>::replicate(Vectors<T, DT, DEVICE> const& vec_in)
{
    this->resize(vec_in.getNumSubs(), vec_in.getShOrder(), 
        vec_in.getGridDim());
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
Scalars<T, DT, DEVICE>::iterator Scalars<T, DT, DEVICE>::getSubN(size_t n)
{
    return(data_ + n * this->getTheDim() * this->getStride());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::const_iterator Scalars<T, DT, DEVICE>::getSubN(
    size_t n) const
{
    return(data_ + n * this->getTheDim() * this->getStride());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::iterator Scalars<T, DT, DEVICE>::end()
{
    return(data_ + this->size());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::const_iterator Scalars<T, DT, DEVICE>::end() const
{
    return(data_ + this->size());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
std::ostream& operator<<(std::ostream& output, const Scalars<T, DT, DEVICE> &sc)
{
    output<<"=====================================================\n"
          <<"SH order            : "<<sc.getShOrder()<<"\n"
          <<"Grid size           : ( "
          <<sc.getGridDim().first<<" , "<<sc.getGridDim().second<<" )"<<"\n"
          <<"stride              : "<<sc.getStride()<<"\n"
          <<"Number of functions : "<<sc.getNumSubs()<<"\n"
          <<"=====================================================\n";
    
    for(typename Scalars<T,DT,DEVICE>::const_iterator it = sc.begin(); 
        it !=sc.end(); ++it)
    {
        output<<*it<<" ";

        if((it-sc.begin() + 1)%sc.getGridDim().second == 0)
            output<<endl;
        
        if((it-sc.begin() + 1)%sc.getStride() == 0)
            output<<endl;
    }

    return(output);
}
