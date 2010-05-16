template<typename T, enum DeviceType DT>
const T SHTrans<T,DT>::alpha_;

template<typename T, enum DeviceType DT>
const T SHTrans<T,DT>::beta_;

template<typename T, enum DeviceType DT>
SHTrans<T,DT>::SHTrans(Device<DT> *dev, int p_in) :
    device_(dev),
    mats_(device_, p_in),
    p(p_in),
    dft_size(2*p)
{}

template<typename T, enum DeviceType DT>
SHTrans<T,DT>::~SHTrans() 
{}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::DLT(T *trans, const T *inputs, T *outputs,
    int m, int n , int k, int mf, int nf, int kf) {
    for (int freq = 0; freq <= p; freq++) {
        int num_legendre_inputs = n;
        if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

        device_->gemm("N", "N", &m, &num_legendre_inputs, &k, &alpha_, 
            trans, &m, inputs, &k, &beta_,outputs, &m);

        trans += m * k;
        inputs += num_legendre_inputs * k;
        outputs += m * num_legendre_inputs;
        if (mf) m--;
        if (nf) n--;
        if (kf) k--;
    }
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::back(const T *inputs, T *work_arr, int n_funs,
                    T *outputs, T *trans, T *dft) {
    T *trans_in = work_arr;
    T *trans_out = work_arr + 2 * p * (p + 1) * n_funs;
    int num_dft_inputs = n_funs * (p + 1);
    DLT(trans, inputs, trans_in, p + 1, 2 * n_funs, p + 1, 0, 0, 1);
    
    device_->Transpose(trans_in, dft_size, num_dft_inputs, trans_out);
    device_->gemm("T", "N", &dft_size, &num_dft_inputs, 
        &dft_size, &alpha_, dft, &dft_size,
        trans_out, &dft_size, &beta_, outputs, &dft_size);
}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::forward(const T *inputs, T *work_arr, int n_funs,
                       T *outputs) {
    T *trans_in = work_arr;
    T *trans_out = work_arr + 2 * p * (p + 1) * n_funs;
    int num_dft_inputs = n_funs * (p + 1);
    device_->gemm("N", "N", &dft_size, &num_dft_inputs, &dft_size, 
        &alpha_, mats_.dft_, &dft_size,inputs, &dft_size, &beta_, 
        trans_in, &dft_size);
    device_->Transpose(trans_in, num_dft_inputs, dft_size, trans_out);
    DLT(mats_.dlt_, trans_out, outputs, p + 1, 2 * n_funs, p + 1, 1, 0, 0);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward(const T *inputs, T *work_arr, 
    int n_funs, T *outputs)
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_du(const T *inputs, T *work_arr, 
    int n_funs, T *outputs) 
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d1_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2u(const T *inputs, T *work_arr, 
    int n_funs, T *outputs) 
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d2_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_dv(const T *inputs, T *work_arr, 
    int n_funs, T *outputs)
{
    back(inputs, work_arr, n_funs, outputs,mats_.dlt_inv_, mats_.dft_inv_d1_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2v(const T *inputs, T *work_arr,
    int n_funs, T *outputs)
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_, mats_.dft_inv_d2_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_duv(const T *inputs, T *work_arr,
    int n_funs, T *outputs)
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d1_, mats_.dft_inv_d1_);
}
