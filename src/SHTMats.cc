template<typename T>
SHTMats<T>::SHTMats(int sh_order, T *data, bool generateMats,
    pair<int, int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? gridDimOf(sh_order_) : grid_dim),
    data_(data),
    dft_size(grid_dim_.second)
{
    assert(data_ != NULL);
    
    size_t dft_size = getDFTLength();
    size_t dlt_size = getDLTLength();

    dft_        = data_;
    dft_inv_    = dft_        + dft_size;
    dft_inv_d1_ = dft_inv_    + dft_size;
    dft_inv_d2_ = dft_inv_d1_ + dft_size;
    dlt_        = dft_inv_d2_ + dft_size;

    dlt_inv_    = dlt_        + dlt_size;
    dlt_inv_d1_ = dlt_inv_    + dlt_size;
    dlt_inv_d2_ = dlt_inv_d1_ + dlt_size;

    if(generateMats)
    {
        gen_dft_forward();
        gen_dft_backward();
        gen_dft_d1backward();
        gen_dft_d2backward();
    }
}

template<typename T>
SHTMats<T>::~SHTMats()
{}

template<typename T>
int SHTMats<T>::getShOrder() const
{
    return(sh_order_);
}

template<typename T>
pair<int, int> SHTMats<T>::getGridDim() const
{
    return(grid_dim_);
}

template<typename T>
size_t SHTMats<T>::getDFTLength() const
{
    return(grid_dim_.second * grid_dim_.second);
}

template<typename T>
size_t SHTMats<T>::getDLTLength() const
{
    return(grid_dim_.first * grid_dim_.first * ( grid_dim_.first + 1) );
}

template<typename T>
size_t SHTMats<T>::getDataLength() const 
{
    return( 4 * getDFTLength() + 4 * getDLTLength());
}

template<typename T>
size_t SHTMats<T>::getDataLength(int sh_order, pair<int, int> grid_dim)
{
    grid_dim = ((grid_dim == EMPTY_GRID) ? gridDimOf(sh_order) : grid_dim);
    
    return( 4 * grid_dim.second * grid_dim.second + 
        4 * grid_dim.first * grid_dim.first * ( grid_dim.first + 1));
}

template<typename T>
void SHTMats<T>::gen_dft_forward() {
    for(int j = 0; j < dft_size; j++)
        dft_[dft_size * j] = 1.0F/dft_size;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / sh_order_) 
                / dft_size * 2;
            
            dft_[2 * i + dft_size * j] = sin(M_PI * i * j / sh_order_) 
                / dft_size * 2;
        }

    for(int j = 0; j < dft_size; j++)
        dft_[dft_size - 1 + dft_size * j] = cos(M_PI * j) / dft_size;
}


template<typename T>
void SHTMats<T>::gen_dft_backward() {
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


template<typename T>
void SHTMats<T>::gen_dft_d1backward() {
    for(int j = 0; j < dft_size; j++)
        dft_inv_d1_[dft_size * j] = 0;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_inv_d1_[2 * i - 1 + dft_size * j] = -i * 
                sin(M_PI * i * j / sh_order_);
            
            dft_inv_d1_[2 * i + dft_size * j] = i * 
                cos(M_PI * i * j / sh_order_);
        }

    for(int j = 0; j < dft_size; j++)
        dft_inv_d1_[dft_size - 1 + dft_size * j] = 0;
}


template<typename T>
void SHTMats<T>::gen_dft_d2backward() {
    for(int j = 0; j < dft_size; j++)
        dft_inv_d2_[dft_size * j] = 0;
  
    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            dft_inv_d2_[2 * i - 1 + dft_size * j] = -i * i * 
                cos(M_PI * i * j / sh_order_);

            dft_inv_d2_[2 * i + dft_size * j] = -i * i * 
                sin(M_PI * i * j / sh_order_);
        }

    for(int j=0; j<dft_size; j++)
        dft_inv_d2_[dft_size-1 + dft_size*j] = - sh_order_*sh_order_*cos(M_PI*j);
}
