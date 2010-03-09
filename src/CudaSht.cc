#include "CudaSht.h"

void cuda_sht::gen_dft_forward() {
  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size * j] = 1.0F/dft_size;

  for(int j = 0; j < dft_size; j++)
    for(int i = 1; i < p; i++) {
      dft_temp[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / p) / dft_size * 2;
      dft_temp[2 * i + dft_size * j] = sin(M_PI * i * j / p) / dft_size * 2;
    }

  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size - 1 + dft_size * j] = cos(M_PI * j) / dft_size;

  cublas_alloc_copy(dft_temp, &dft_forward, dft_size, dft_size);
}


void cuda_sht::gen_dft_backward() {
  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size * j] = 1.0F;
  
  for(int j = 0; j < dft_size; j++)
    for(int i = 1; i < p; i++) {
      dft_temp[2 * i - 1 + dft_size * j] = cos(M_PI * i * j / p);
      dft_temp[2 * i + dft_size * j] = sin(M_PI * i * j / p);
    }

  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size - 1 + dft_size * j] = cos(M_PI * j);

  cublas_alloc_copy(dft_temp, &dft_backward, dft_size, dft_size);
}


void cuda_sht::gen_dft_d1backward() {
  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size * j] = 0;

  for(int j = 0; j < dft_size; j++)
    for(int i = 1; i < p; i++) {
      dft_temp[2 * i - 1 + dft_size * j] = -i * sin(M_PI * i * j / p);
      dft_temp[2 * i + dft_size * j] = i * cos(M_PI * i * j / p);
    }

  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size - 1 + dft_size * j] = 0;

  cublas_alloc_copy(dft_temp, &dft_d1backward, dft_size, dft_size);
}


void cuda_sht::gen_dft_d2backward() {
  for(int j = 0; j < dft_size; j++)
    dft_temp[dft_size * j] = 0;
  
  for(int j = 0; j < dft_size; j++)
    for(int i = 1; i < p; i++) {
      dft_temp[2 * i - 1 + dft_size * j] = -i * i * cos(M_PI * i * j / p);
      dft_temp[2 * i + dft_size * j] = -i * i * sin(M_PI * i * j / p);
    }

  for(int j=0; j<dft_size; j++)
    dft_temp[dft_size-1 + dft_size*j] = -p*p*cos(M_PI*j);

  cublas_alloc_copy(dft_temp, &dft_d2backward, dft_size, dft_size);
}


void cuda_sht::cublas_alloc_copy(scalar *cpu_ptr, scalar **gpu_ptr, int rows, int cols) {
  cublasAlloc(rows * cols, sizeof(scalar), (void**) gpu_ptr);
  cublasSetMatrix(rows, cols, sizeof(scalar), (void*) cpu_ptr, rows, (void*) *gpu_ptr, rows);
}


void cuda_sht::read_leg_mat(scalar **gpu_ptr, char *fname) {
  std::ifstream file(fname);
  if (file.is_open()) {
    int idx=0;
    while (idx < leg_mat_size) {
      file >> leg_temp[idx++];
    }
    file.close();
    cublas_alloc_copy(leg_temp, gpu_ptr, (p + 1), (p + 1) * (p + 2));
  } else {
    std::cout << "Unable to open file"; 
    abort();
  }
}


cuda_sht::cuda_sht(int p, int num_vesicles, char *leg_trans_fname, char *leg_trans_inv_fname,
                   char *d1_leg_trans_fname, char *d2_leg_trans_fname, scalar *leg_temp,
                   scalar *dft_temp, scalar *trans_in_cpu, scalar *trans_out_cpu) {
  this->p = p;
  this->num_vesicles = num_vesicles;
  this->dft_size = 2 * p;
  this->num_dft_inputs = num_vesicles * (p + 1);
  this->leg_mat_size = (p + 1) * (p + 1) * (p + 2);
  this->vesicle_size = 2 * p * (p + 1);
  this->leg_temp = leg_temp;
  this->trans_in = trans_in;
  this->trans_out = trans_out;
  this->dft_temp = dft_temp;
  this->trans_in_cpu = trans_in_cpu;
  this->trans_out_cpu = trans_out_cpu;

  cublasInit();
  cublasAlloc(2 * p * (p + 1) * num_vesicles, sizeof(scalar), (void**) &trans_in);
  cublasAlloc(2 * p * (p + 1) * num_vesicles, sizeof(scalar), (void**) &trans_out);

  read_leg_mat(&leg_trans, leg_trans_fname);
  read_leg_mat(&leg_trans_inv, leg_trans_inv_fname);
  read_leg_mat(&d1_leg_trans, d1_leg_trans_fname);
  read_leg_mat(&d2_leg_trans, d2_leg_trans_fname);
  gen_dft_forward();
  gen_dft_backward();
  gen_dft_d1backward();
  gen_dft_d2backward();
}


cuda_sht::~cuda_sht() {
  cublasShutdown();
}


void cuda_sht::pre_legendre(scalar *in, scalar *out) {
  cublasGetMatrix(dft_size, num_dft_inputs, sizeof(scalar), in, dft_size, trans_in_cpu, dft_size);

  int leg_input_pointer = 0;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++) {
        trans_out_cpu[leg_input_pointer++] =
            trans_in_cpu[vesicle_size*vesicle_num + dft_size*i + freq];
      }

  cublasSetMatrix(dft_size, num_dft_inputs, sizeof(scalar), (void*) trans_out_cpu,
                  dft_size, (void*) out, dft_size);
}


void cuda_sht::post_legendre(scalar *in, scalar *out) {
  cublasGetMatrix(dft_size, num_dft_inputs, sizeof(scalar), in, dft_size, trans_in_cpu, dft_size);

  int leg_input_pointer = 0;
  for (int freq=0; freq<dft_size; freq++)
    for (int vesicle_num=0; vesicle_num<num_vesicles; vesicle_num++)
      for (int i=0; i<p+1; i++)
        trans_out_cpu[vesicle_size * vesicle_num + dft_size * i + freq] = trans_in_cpu[leg_input_pointer++];

  cublasSetMatrix(dft_size, num_dft_inputs, sizeof(scalar), (void*) trans_out_cpu,
                  dft_size, (void*) out, dft_size);
}


void cuda_sht::leg_transform(scalar *trans_gpu, scalar *inputs_gpu, scalar *outputs_gpu,
                int m, int n , int k, int mf, int nf, int kf) {
  for (int freq = 0; freq <= p; freq++) {
    int num_legendre_inputs = n;
    if (freq == 0 || freq == p) num_legendre_inputs = n / 2;

    cublasSgemm('N', 'N', m, num_legendre_inputs, k, 1.0F, trans_gpu, m, inputs_gpu, k, 0.0F,
                outputs_gpu, m);

    trans_gpu += m * k;
    inputs_gpu += num_legendre_inputs * k;
    outputs_gpu += m * num_legendre_inputs;
    if (mf) m--;
    if (nf) n--;
    if (kf) k--;
  }
}


void cuda_sht::back(scalar *inputs, scalar *outputs, scalar *trans, scalar *dft) {
  leg_transform(trans, inputs, trans_in, p + 1, 2 * num_vesicles, p + 1, 0, 0, 1);
  cu_trans(trans_out, trans_in, num_dft_inputs, dft_size);
  cublasSgemm('T', 'N', dft_size, num_dft_inputs, dft_size, 1.0F, dft,dft_size,
              trans_out, dft_size, 0.0F, outputs, dft_size);
}

void cuda_sht::forward(scalar *inputs, scalar *outputs) {
  cublasSgemm('N', 'N', dft_size, num_dft_inputs, dft_size, 1.0F, dft_forward, dft_size,
              inputs, dft_size, 0.0F, trans_in, dft_size);
  cu_trans(trans_out, trans_in, dft_size, num_dft_inputs);
  leg_transform(leg_trans, trans_out, outputs, p + 1, 2 * num_vesicles, p + 1, 1, 0, 0);
}


void cuda_sht::backward(scalar *inputs, scalar *outputs) {
  back(inputs, outputs, leg_trans_inv, dft_backward);
}


void cuda_sht::backward_du(scalar *inputs, scalar *outputs) {
  back(inputs, outputs, d1_leg_trans, dft_backward);
}


void cuda_sht::backward_d2u(scalar *inputs, scalar *outputs) {
  back(inputs, outputs, d2_leg_trans, dft_backward);
}


void cuda_sht::backward_dv(scalar *inputs, scalar *outputs) {
  back(inputs, outputs,leg_trans_inv, dft_d1backward);
}


void cuda_sht::backward_d2v(scalar *inputs, scalar *outputs) {
  back(inputs, outputs, leg_trans_inv, dft_d2backward);
}


void cuda_sht::backward_duv(scalar *inputs, scalar *outputs) {
  back(inputs, outputs, d1_leg_trans, dft_d1backward);
}


void cuda_sht::test() {
  scalar *inputs, *outputs, *outputs_2;
  cublasAlloc(2 * p * (p + 1) * num_vesicles, sizeof(scalar), (void**) &outputs);
  cublasAlloc(2 * p * (p + 1) * num_vesicles, sizeof(scalar), (void**) &outputs_2);
  float* input_cpu = (float*)malloc(sizeof(float) * 2 * p * (p + 1) * num_vesicles);
  float* output_cpu = (float*)malloc(sizeof(float) * 2 * p * (p + 1) * num_vesicles);
  float* output_cpu_2 = (float*)calloc(sizeof(float), 2 * p * (p + 1) * num_vesicles);

  for (int i = 0; i < 2 * p * (p + 1) * num_vesicles; i++) {
    input_cpu[i] = (float) rand() / (float) RAND_MAX;
  }
  cublas_alloc_copy(input_cpu, &inputs, dft_size, num_dft_inputs);
  unsigned int timer;
  cutCreateTimer(&timer);
  cutStartTimer(timer);
  forward(inputs, outputs);
  backward(outputs, outputs_2);
  cutStopTimer(timer);
  fprintf(stderr, "forward backward took %fms\n", cutGetTimerValue(timer));
  cublasGetMatrix(dft_size, num_dft_inputs, sizeof(scalar), outputs_2, dft_size,
                  output_cpu_2, dft_size);
  for (int i = 0; i < 2* p * (p + 1) * num_vesicles; i++) {
    fprintf(stderr, "%f\n", output_cpu_2[i]);
  }
}
