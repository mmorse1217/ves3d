#include "BlasSht.h"

void BlasSht::gen_dft_forward() {
    for(int j = 0; j < dft_size; j++)
        dft_forward[dft_size * j] = 1.0F/dft_size;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < p; i++) {
            dft_forward[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / p) / dft_size * 2;
            dft_forward[2 * i + dft_size * j] = sin(M_PI * i * j / p) / dft_size * 2;
        }

    for(int j = 0; j < dft_size; j++)
        dft_forward[dft_size - 1 + dft_size * j] = cos(M_PI * j) / dft_size;
}


void BlasSht::gen_dft_backward() {
    for(int j = 0; j < dft_size; j++)
        dft_backward[dft_size * j] = 1.0F;
  
    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < p; i++) {
            dft_backward[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / p);
            dft_backward[2 * i + dft_size * j] = sin(M_PI * i * j / p);
        }

    for(int j = 0; j < dft_size; j++)
        dft_backward[dft_size - 1 + dft_size * j] = cos(M_PI * j);
}


void BlasSht::gen_dft_d1backward() {
    for(int j = 0; j < dft_size; j++)
        dft_d1backward[dft_size * j] = 0;

    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < p; i++) {
            dft_d1backward[2 * i - 1 + dft_size * j] = -i * sin(M_PI * i * j / p);
            dft_d1backward[2 * i + dft_size * j] = i * cos(M_PI * i * j / p);
        }

    for(int j = 0; j < dft_size; j++)
        dft_d1backward[dft_size - 1 + dft_size * j] = 0;
}


void BlasSht::gen_dft_d2backward() {
    for(int j = 0; j < dft_size; j++)
        dft_d2backward[dft_size * j] = 0;
  
    for(int j = 0; j < dft_size; j++)
        for(int i = 1; i < p; i++) {
            dft_d2backward[2 * i - 1 + dft_size * j] = -i * i * cos(M_PI * i * j / p);
            dft_d2backward[2 * i + dft_size * j] = -i * i * sin(M_PI * i * j / p);
        }

    for(int j=0; j<dft_size; j++)
        dft_d2backward[dft_size-1 + dft_size*j] = -p*p*cos(M_PI*j);
}


void BlasSht::read_leg_mat(scalar *leg_ptr, char *fname) {
    std::ifstream file(fname);
    if (file.is_open()) {
        int idx=0;
        while (idx < leg_mat_size) {
            file >> leg_ptr[idx++];
        }
        file.close();
    } else {
        std::cout << "Unable to open file"; 
        abort();
    }
}


void BlasSht::transpose(scalar *out, scalar *in, int width, int height) {
    int leg_input_pointer = 0;
    for (int freq=0; freq<width; freq++)
        for (int i=0; i<height; i++) {
            out[leg_input_pointer++] = in[i * width + freq];
        }
}


BlasSht::BlasSht() :
    p(0),
    dft_size(0),
    leg_mat_size(0),
    vesicle_size(0),
    alpha(0),
    beta(0),
    trans_in(0),
    trans_out(0),
    leg_trans(0),
    leg_trans_inv(0),
    d1_leg_trans(0),
    d2_leg_trans(0),
    dft_forward(0),
    dft_backward(0),
    dft_d1backward(0),
    dft_d2backward(0)
{}

BlasSht::BlasSht(int p, char *leg_trans_fname, char *leg_trans_inv_fname,
    char *d1_leg_trans_fname, char *d2_leg_trans_fname,
    scalar *dft_forward, scalar *dft_backward,
    scalar *dft_d1backward, scalar *dft_d2backward,
    scalar *leg_trans, scalar *leg_trans_inv,
    scalar* d1_leg_trans, scalar *d2_leg_trans,
    scalar *trans_in, scalar *trans_out) {
    
    InitializeBlasSht(p, leg_trans_fname, leg_trans_inv_fname, 
        d1_leg_trans_fname, d2_leg_trans_fname, dft_forward, dft_backward, dft_d1backward, 
        dft_d2backward, leg_trans, leg_trans_inv, d1_leg_trans, d2_leg_trans, trans_in, trans_out);

//     this->p = p;
//     //this->num_vesicles = num_vesicles;
//     this->dft_size = 2 * p;
//     //this->num_dft_inputs = num_vesicles * (p + 1);
//     this->leg_mat_size = (p + 1) * (p + 1) * (p + 2);
//     this->vesicle_size = 2 * p * (p + 1);
//     this->trans_in = trans_in;
//     this->trans_out = trans_out;
//     this->dft_forward = dft_forward;
//     this->dft_backward = dft_backward;
//     this->dft_d1backward = dft_d1backward;
//     this->dft_d2backward = dft_d2backward;
//     this->leg_trans = leg_trans;
//     this->leg_trans_inv = leg_trans_inv;
//     this->d1_leg_trans = d1_leg_trans;
//     this->d2_leg_trans= d2_leg_trans;
//     this->alpha = 1.0F;
//     this->beta = 0.0F;

//     read_leg_mat(leg_trans, leg_trans_fname);
//     read_leg_mat(leg_trans_inv, leg_trans_inv_fname);
//     read_leg_mat(d1_leg_trans, d1_leg_trans_fname);
//     read_leg_mat(d2_leg_trans, d2_leg_trans_fname);
//     gen_dft_forward();
//     gen_dft_backward();
//     gen_dft_d1backward();
//     gen_dft_d2backward();
}


BlasSht::~BlasSht() {
}


void BlasSht::leg_transform(scalar *trans, const scalar *inputs, scalar *outputs,
    int m, int n , int k, int mf, int nf, int kf) {
    for (int freq = 0; freq <= p; freq++) {
        int num_legendre_inputs = n;
        if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

        sgemm("N", "N", &m, &num_legendre_inputs, &k, &alpha, trans, &m, inputs, &k, &beta,
            outputs, &m);

        trans += m * k;
        inputs += num_legendre_inputs * k;
        outputs += m * num_legendre_inputs;
        if (mf) m--;
        if (nf) n--;
        if (kf) k--;
    }
}


void BlasSht::back(const scalar *inputs, int n_funs, scalar *outputs, scalar *trans, scalar *dft) {

    int num_dft_inputs = n_funs * (p + 1);
    leg_transform(trans, inputs, trans_in, p + 1, 2 * n_funs, p + 1, 0, 0, 1);
    transpose(trans_out, trans_in, num_dft_inputs, dft_size);
    sgemm("T", "N", &dft_size, &num_dft_inputs, &dft_size, &alpha, dft, &dft_size,
        trans_out, &dft_size, &beta, outputs, &dft_size);
}

void BlasSht::forward(const scalar *inputs, int n_funs, scalar *outputs) {

    int num_dft_inputs = n_funs * (p + 1);
    sgemm("N", "N", &dft_size, &num_dft_inputs, &dft_size, &alpha, dft_forward, &dft_size,
        inputs, &dft_size, &beta, trans_in, &dft_size);
    transpose(trans_out, trans_in, dft_size, num_dft_inputs);
    leg_transform(leg_trans, trans_out, outputs, p + 1, 2 * n_funs, p + 1, 1, 0, 0);
}


void BlasSht::backward(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs, leg_trans_inv, dft_backward);
}


void BlasSht::backward_du(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs, d1_leg_trans, dft_backward);
}


void BlasSht::backward_d2u(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs, d2_leg_trans, dft_backward);
}


void BlasSht::backward_dv(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs,leg_trans_inv, dft_d1backward);
}


void BlasSht::backward_d2v(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs, leg_trans_inv, dft_d2backward);
}


void BlasSht::backward_duv(const scalar *inputs, int n_funs, scalar *outputs) {
    back(inputs, n_funs, outputs, d1_leg_trans, dft_d1backward);
}


void BlasSht::test(int n_funs) {
    scalar *inputs, *outputs, *outputs_2;
    inputs = (scalar*) malloc(2 * p * (p + 1) * n_funs * sizeof(scalar));
    outputs = (scalar*) malloc(2 * p * (p + 1) * n_funs * sizeof(scalar));
    outputs_2 = (scalar*) malloc(2 * p * (p + 1) * n_funs * sizeof(scalar));

    for (int i = 0; i < 2 * p * (p + 1) * n_funs; i++) {
        inputs[i] = (float) rand() / (float) RAND_MAX;
    }
    forward(inputs, n_funs, outputs);
    backward(outputs, n_funs, outputs_2);
//     for (int i = 0; i < 2* p * (p + 1); i++) {
//         fprintf(stderr, "%f\n", outputs_2[i]);
//     }
}

void BlasSht::InitializeBlasSht(int p, char *leg_trans_fname, 
    char *leg_trans_inv_fname, char *d1_leg_trans_fname, char *d2_leg_trans_fname,
    scalar *dft_forward, scalar *dft_backward, scalar *dft_d1backward, 
    scalar *dft_d2backward, scalar *leg_trans, scalar *leg_trans_inv, 
    scalar* d1_leg_trans, scalar *d2_leg_trans, scalar *trans_in, scalar *trans_out) 
{
    this->p = p;
    this->dft_size = 2 * p;
    //this->num_dft_inputs = num_vesicles * (p + 1);
    this->leg_mat_size = (p + 1) * (p + 1) * (p + 2);
    this->vesicle_size = 2 * p * (p + 1);
    this->trans_in = trans_in;
    this->trans_out = trans_out;
    this->dft_forward = dft_forward;
    this->dft_backward = dft_backward;
    this->dft_d1backward = dft_d1backward;
    this->dft_d2backward = dft_d2backward;
    this->leg_trans = leg_trans;
    this->leg_trans_inv = leg_trans_inv;
    this->d1_leg_trans = d1_leg_trans;
    this->d2_leg_trans= d2_leg_trans;
    this->alpha = 1.0F;
    this->beta = 0.0F;

    read_leg_mat(leg_trans, leg_trans_fname);
    read_leg_mat(leg_trans_inv, leg_trans_inv_fname);
    read_leg_mat(d1_leg_trans, d1_leg_trans_fname);
    read_leg_mat(d2_leg_trans, d2_leg_trans_fname);
    
    gen_dft_forward();
    gen_dft_backward();
    gen_dft_d1backward();
    gen_dft_d2backward();
}

