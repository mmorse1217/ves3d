template<typename T, enum DeviceType DT>
const T SHTrans<T,DT>::alpha_;

template<typename T, enum DeviceType DT>
const T SHTrans<T,DT>::beta_;

template<typename T, enum DeviceType DT>
SHTrans<T,DT>::SHTrans(const Device<DT> *dev, int p_in) :
    device_(dev),
    mats_(device_, p_in),
    p(p_in),
    dft_size(2*p),
    filter_coeff_((T*) device_->Malloc(p * (p + 2) * sizeof(T)))
{
    int ll = p * (p + 2);
    T *buffer = (T*) malloc(p * (p + 2) * sizeof(T));
    
    int idx = 0, len;
    for(int ii=0; ii< 2 * p; ++ii)
    {
        len = p + 1 - (ii+1)/2;
        for(int jj=0; jj < len; ++jj)
             buffer[idx++] = (len-jj)<=(p-2*p/3) ? 0 : 1;
    }
    device_->Memcpy(filter_coeff_, buffer, p *(p + 2) * sizeof(T), 
        MemcpyHostToDevice);
    free(buffer);
}

template<typename T, enum DeviceType DT>
SHTrans<T,DT>::~SHTrans() 
{
    device_->Free(filter_coeff_);
}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::DLT(T *trans, const T *inputs, T *outputs,
    int m, int n , int k, int mf, int nf, int kf) const 
{
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
                    T *outputs, T *trans, T *dft) const 
{
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
                       T *outputs) const
{
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
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_du(const T *inputs, T *work_arr, 
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d1_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2u(const T *inputs, T *work_arr, 
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d2_, mats_.dft_inv_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_dv(const T *inputs, T *work_arr, 
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs,mats_.dlt_inv_, mats_.dft_inv_d1_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2v(const T *inputs, T *work_arr,
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_, mats_.dft_inv_d2_);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_duv(const T *inputs, T *work_arr,
    int n_funs, T *outputs) const
{
    back(inputs, work_arr, n_funs, outputs, mats_.dlt_inv_d1_, mats_.dft_inv_d1_);
}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::Filter(const T *inputs, T *work_arr, int n_funs, T* shc, T *outputs) const
{
    forward(inputs, work_arr, n_funs, shc);
    ScaleFreq(shc, n_funs, this->filter_coeff_, shc);
    backward(shc, work_arr, n_funs, outputs);
}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::ScaleFreq(const T *shc_in, int n_funs, const T* scaling_coeff, T *shc_out) const
{
    int leg_order = p+1;
    
    device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
    scaling_coeff += leg_order;
    shc_in += n_funs * leg_order;
    shc_out += n_funs * leg_order;
    leg_order--;

    // process remaining frequencies except the last cosine
    for (; leg_order>1; leg_order--) 
    {
        // first process cosine
        device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
        
        // then process sine
        device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
        scaling_coeff += leg_order;
        shc_in += n_funs * leg_order;
        shc_out += n_funs * leg_order;
    }
    
    // process last cosine
    device_->ax(scaling_coeff, shc_in, leg_order, n_funs, shc_out);
}
