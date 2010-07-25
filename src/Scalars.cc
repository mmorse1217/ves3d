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
    num_subs_(num_subs)
{
    Array<T, DT, DEVICE>::resize(this->getStride() * this->getTheDim() 
        * this->getNumSubs());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
Scalars<T, DT, DEVICE>::~Scalars()
{}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
int Scalars<T, DT, DEVICE>::getTheDim()
{
    return(the_dim_);
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
size_t Scalars<T, DT, DEVICE>::getNumSubs() const
{
    return(num_subs_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
size_t Scalars<T, DT, DEVICE>::getSubLength() const
{
    return((num_subs_ > 0 ) * the_dim_ * stride_);
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
void Scalars<T, DT, DEVICE>::resize(size_t new_num_subs, int new_sh_order,
    pair<int, int> new_grid_dim)
{
    num_subs_ = new_num_subs;
    sh_order_ = (new_sh_order == -1) ? sh_order_ : new_sh_order;
    grid_dim_ = (new_grid_dim == EMPTY_GRID) ? 
        gridDimOf(sh_order_) : new_grid_dim;
    stride_   = grid_dim_.first * grid_dim_.second;

    Array<T, DT, DEVICE>::resize(this->getStride() * this->getTheDim() 
        * this->getNumSubs());

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
T* Scalars<T, DT, DEVICE>::getSubN(size_t n)
{
    return(Array<T, DT, DEVICE>::begin() + n * this->getTheDim()
        * this->getStride());
}

template<typename T, enum DeviceType DT, const Device<DT> &DEVICE> 
const T* Scalars<T, DT, DEVICE>::getSubN(
    size_t n) const
{
    return(Array<T, DT, DEVICE>::begin() + n * this->getTheDim() 
        * this->getStride());
}

/////////////////////////////////////////////////////////////////////////////////
template<typename T, enum DeviceType DT, const Device<DT> &DEVICE>
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
    this->c_.getDevice().Memcpy(buffer, this->c_.begin(), 
        this->c_.size() * sizeof(T), MemcpyDeviceToHost);
    
    for(size_t ii=0; ii<this->c_.size(); ++ii)
    {
        output<<buffer[ii]<<" ";
        
        if((ii + 1)%this->c_.getGridDim().second == 0)
            output<<endl;
        
        if((ii + 1)%this->c_.getStride() == 0)
            output<<endl;
    }

    delete[] buffer;
    
    return(output);
}

template<typename Container>
std::ostream& operator<<(std::ostream& output, const ShowEntries<Container> &se)
{
    return(se(output));
}
    
