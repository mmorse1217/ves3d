template<typename T, enum DeviceType DT>
SHTMats<T, DT>::SHTMats(const Device<DT> *dev, int sh_order,
    OperatorsMats<T> &mats, pair<int,int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? gridDimOf(sh_order_) : grid_dim),
    device_(dev),
    data_((T*) device_->Malloc(this->GetDataLength() * sizeof(T))),
    dft_size(grid_dim_.second)
{
    size_t dft_size = GetDFTLength();
    size_t dlt_size = GetDLTLength();

    dft_        = data_;
    dft_inv_    = dft_        + dft_size;
    dft_inv_d1_ = dft_inv_    + dft_size;
    dft_inv_d2_ = dft_inv_d1_ + dft_size;
    dlt_        = dft_inv_d2_ + dft_size;

    dlt_inv_    = dlt_        + dlt_size;
    dlt_inv_d1_ = dlt_inv_    + dlt_size;
    dlt_inv_d2_ = dlt_inv_d1_ + dlt_size;

    gen_dft_forward();
    gen_dft_backward();
    gen_dft_d1backward();
    gen_dft_d2backward();

    ///@todo this should be computed rather than read from file
    DataIO<T,DT> IO(*device_,"",0);
    char fname[500];

    
    sprintf(fname,"precomputed/legTrans%u_single.txt",sh_order_);
    IO.ReadData(fname, dlt_size, dlt_);
    
    sprintf(fname,"precomputed/legTransInv%u_single.txt",sh_order_);
    IO.ReadData(fname, dlt_size, dlt_inv_);
    
    sprintf(fname,"precomputed/d1legTrans%u_single.txt",sh_order_);
    IO.ReadData(fname, dlt_size, dlt_inv_d1_);
    
    sprintf(fname,"precomputed/d2legTrans%u_single.txt",sh_order_);
    IO.ReadData(fname, dlt_size, dlt_inv_d2_);
}

template<typename T, enum DeviceType DT>
SHTMats<T, DT>::~SHTMats()
{
    device_->Free(data_);
}

template<typename T, enum DeviceType DT>
int SHTMats<T, DT>::GetShOrder() const
{
    return(sh_order_);
}

template<typename T, enum DeviceType DT>
pair<int, int> SHTMats<T, DT>::GetGridDim() const
{
    return(grid_dim_);
}

template<typename T, enum DeviceType DT>
const Device<DT>* SHTMats<T, DT>::GetDevicePtr() const
{
    return(device_);
}

template<typename T, enum DeviceType DT>
size_t SHTMats<T, DT>::GetDFTLength() const
{
    return(grid_dim_.second * grid_dim_.second);
}

template<typename T, enum DeviceType DT>
size_t SHTMats<T, DT>::GetDLTLength() const
{
    return(grid_dim_.first * grid_dim_.first * ( grid_dim_.first + 1) );
}

template<typename T, enum DeviceType DT>
size_t SHTMats<T, DT>::GetDataLength() const
{
    return( 4 * GetDFTLength() + 4 * GetDLTLength());
}

template<typename T, enum DeviceType DT>
void SHTMats<T,DT>::gen_dft_forward() {
    for(int j = 0; j < dft_size; j++)
        dft_[dft_size * j] = 1.0F/dft_size;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / sh_order_) / dft_size * 2;
            dft_[2 * i + dft_size * j] = sin(M_PI * i * j / sh_order_) / dft_size * 2;
        }

    for(int j = 0; j < dft_size; j++)
        dft_[dft_size - 1 + dft_size * j] = cos(M_PI * j) / dft_size;
}


template<typename T, enum DeviceType DT>
void SHTMats<T,DT>::gen_dft_backward() {
    for(int j = 0; j < dft_size; j++)
        dft_inv_[dft_size * j] = 1.0F;
  
    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_inv_[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / sh_order_);
            dft_inv_[2 * i + dft_size * j] = sin(M_PI * i * j / sh_order_);
        }

    for(int j = 0; j < dft_size; j++)
        dft_inv_[dft_size - 1 + dft_size * j] = cos(M_PI * j);
}


template<typename T, enum DeviceType DT>
void SHTMats<T,DT>::gen_dft_d1backward() {
    for(int j = 0; j < dft_size; j++)
        dft_inv_d1_[dft_size * j] = 0;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_inv_d1_[2 * i - 1 + dft_size * j] = -i * sin(M_PI * i * j / sh_order_);
            dft_inv_d1_[2 * i + dft_size * j] = i * cos(M_PI * i * j / sh_order_);
        }

    for(int j = 0; j < dft_size; j++)
        dft_inv_d1_[dft_size - 1 + dft_size * j] = 0;
}


template<typename T, enum DeviceType DT>
void SHTMats<T,DT>::gen_dft_d2backward() {
    for(int j = 0; j < dft_size; j++)
        dft_inv_d2_[dft_size * j] = 0;
  
    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_inv_d2_[2 * i - 1 + dft_size * j] = -i * i * cos(M_PI * i * j / sh_order_);
            dft_inv_d2_[2 * i + dft_size * j] = -i * i * sin(M_PI * i * j / sh_order_);
        }

    for(int j=0; j<dft_size; j++)
        dft_inv_d2_[dft_size-1 + dft_size*j] = - sh_order_*sh_order_*cos(M_PI*j);
}
