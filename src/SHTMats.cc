template<typename T, typename Device>
SHTMats<T, Device>::SHTMats(const Device &dev, int sh_order, T *data,
    bool generateMats, std::pair<int, int> grid_dim) :
    sh_order_(sh_order),
    grid_dim_((grid_dim == EMPTY_GRID) ? SpharmGridDim(sh_order_) : grid_dim),
    data_(data),
    dft_size(grid_dim_.second),
    device_(dev)
{
    ASSERT(data_ != NULL,"NULL pointer passed!");

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

template<typename T, typename Device>
SHTMats<T, Device>::~SHTMats()
{}

template<typename T, typename Device>
int SHTMats<T, Device>::getShOrder() const
{
    return(sh_order_);
}

template<typename T, typename Device>
std::pair<int, int> SHTMats<T, Device>::getGridDim() const
{
    return(grid_dim_);
}

template<typename T, typename Device>
size_t SHTMats<T, Device>::getDFTLength() const
{
    return(grid_dim_.second * grid_dim_.second);
}

template<typename T, typename Device>
size_t SHTMats<T, Device>::getDLTLength() const
{
    return(grid_dim_.first * grid_dim_.first * ( grid_dim_.first + 1) );
}

template<typename T, typename Device>
size_t SHTMats<T, Device>::getDataLength() const
{
    return( 4 * getDFTLength() + 4 * getDLTLength());
}

template<typename T, typename Device>
size_t SHTMats<T, Device>::getDataLength(int sh_order, std::pair<int, int> grid_dim)
{
    grid_dim = ((grid_dim == EMPTY_GRID) ? SpharmGridDim(sh_order) : grid_dim);

    return( 4 * grid_dim.second * grid_dim.second +
        4 * grid_dim.first * grid_dim.first * ( grid_dim.first + 1));
}

template<typename T, typename Device>
const Device& SHTMats<T,Device>::getDevice() const
{
    return(device_);
}

template<typename T, typename Device>
void SHTMats<T, Device>::gen_dft_forward() {

    T* buffer = new T[getDFTLength()];

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size * j] = 1.0F/dft_size;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            buffer[2 * i - 1 + dft_size * j] =
                cos(PI64<T>() * i * j / sh_order_)/ dft_size * 2;

            buffer[2 * i + dft_size * j] =
                sin(PI64<T>() * i * j / sh_order_)/ dft_size * 2;
        }

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size - 1 + dft_size * j] =
            cos(PI64<T>() * j) / dft_size;

    device_.Memcpy(dft_, buffer, getDFTLength() * sizeof(T),
        Device::MemcpyHostToDevice);
    delete[] buffer;
}


template<typename T, typename Device>
void SHTMats<T, Device>::gen_dft_backward() {

    T* buffer = new T[getDFTLength()];

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size * j] = 1.0F;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            buffer[2 * i - 1 + dft_size * j] =
                cos(PI64<T>() * i * j / sh_order_);
            buffer[2 * i + dft_size * j] =
                sin(PI64<T>() * i * j / sh_order_);
        }

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size - 1 + dft_size * j] = cos(PI64<T>() * j);

    device_.Memcpy(dft_inv_, buffer, getDFTLength() * sizeof(T),
        Device::MemcpyHostToDevice);
    delete[] buffer;
}


template<typename T, typename Device>
void SHTMats<T, Device>::gen_dft_d1backward() {

    T* buffer = new T[getDFTLength()];

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size * j] = 0;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            buffer[2 * i - 1 + dft_size * j] = -i *
                sin(PI64<T>() * i * j / sh_order_);

            buffer[2 * i + dft_size * j] = i *
                cos(PI64<T>() * i * j / sh_order_);
        }

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size - 1 + dft_size * j] = 0;

    device_.Memcpy(dft_inv_d1_, buffer, getDFTLength() * sizeof(T),
        Device::MemcpyHostToDevice);
    delete[] buffer;
}


template<typename T, typename Device>
void SHTMats<T, Device>::gen_dft_d2backward() {

    T* buffer = new T[getDFTLength()];

    for(int j = 0; j < dft_size; j++)
        buffer[dft_size * j] = 0;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < sh_order_; i++) {
            buffer[2 * i - 1 + dft_size * j] = -i * i *
                cos(PI64<T>() * i * j / sh_order_);

            buffer[2 * i + dft_size * j] = -i * i *
                sin(PI64<T>() * i * j / sh_order_);
        }

    for(int j=0; j<dft_size; j++)
        buffer[dft_size-1 + dft_size*j] = -
            sh_order_*sh_order_*cos(PI64<T>()*j);

    device_.Memcpy(dft_inv_d2_,
        buffer,
        getDFTLength() * sizeof(T),
        Device::MemcpyHostToDevice);
    delete[] buffer;
}
