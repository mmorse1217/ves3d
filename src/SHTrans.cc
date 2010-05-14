template<typename T, enum DeviceType DT>
SHTrans<T,DT>::SHTrans(Device<DT> *dev, int p_in) :
    device_(dev),
    mats_(dev, p_in),
    p(p_in),
    dft_size(2*p),
    leg_mat_size((p + 1) * (p + 1) * (p + 2)),
    alpha_(1.0),
    beta_(0.0),
    leg_trans(0),
    leg_trans_inv(0),
    d1_leg_trans(0),
    d2_leg_trans(0),
    dft_forward(0),
    dft_backward(0),
    dft_d1backward(0),
    dft_d2backward(0)
{
    this->dft_forward    = mats_.dft_;
    this->dft_backward   = mats_.dft_inv_;
    this->dft_d1backward = mats_.dft_inv_d1_; 
    this->dft_d2backward = mats_.dft_inv_d2_;
                               
    this->leg_trans      = mats_.dlt_;          
    this->leg_trans_inv  = mats_.dlt_inv_;    
    this->d1_leg_trans   = mats_.dlt_inv_d1_;
    this->d2_leg_trans   = mats_.dlt_inv_d2_; 
}

template<typename T, enum DeviceType DT>
SHTrans<T,DT>::~SHTrans() 
{}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::DLT(T *trans, const T *inputs, T *outputs,
    int m, int n , int k, int mf, int nf, int kf) {
    for (int freq = 0; freq <= p; freq++) {
        int num_legendre_inputs = n;
        if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

        device_->gemm("N", "N", &m, &num_legendre_inputs, &k, &alpha_, trans, &m, inputs, &k, &beta_,
            outputs, &m);

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
    device_->gemm("T", "N", &dft_size, &num_dft_inputs, &dft_size, &alpha_, dft, &dft_size,
        trans_out, &dft_size, &beta_, outputs, &dft_size);
}

template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::forward(const T *inputs, T *work_arr, int n_funs,
                       T *outputs) {
    T *trans_in = work_arr;
    T *trans_out = work_arr + 2 * p * (p + 1) * n_funs;
    int num_dft_inputs = n_funs * (p + 1);
    device_->gemm("N", "N", &dft_size, &num_dft_inputs, &dft_size, &alpha_, dft_forward, &dft_size,
        inputs, &dft_size, &beta_, trans_in, &dft_size);
    device_->Transpose(trans_in, num_dft_inputs, dft_size, trans_out);
    DLT(leg_trans, trans_out, outputs, p + 1, 2 * n_funs, p + 1, 1, 0, 0);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs, leg_trans_inv, dft_backward);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_du(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs, d1_leg_trans, dft_backward);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2u(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs, d2_leg_trans, dft_backward);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_dv(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs,leg_trans_inv, dft_d1backward);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_d2v(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs, leg_trans_inv, dft_d2backward);
}


template<typename T, enum DeviceType DT>
void SHTrans<T,DT>::backward_duv(const T *inputs, T *work_arr, int n_funs, T *outputs) {
    back(inputs, work_arr, n_funs, outputs, d1_leg_trans, dft_d1backward);
}
