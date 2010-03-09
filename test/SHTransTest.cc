#include<iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SHTrans.h"


int main(int argc, char ** argv)
{
    const int p = 12;
    const int num_vesicles = 300;
    int leg_mat_size      = (p + 1) * (p + 1) * (p + 2);
    float* dft_forward    = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_backward   = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_d1backward = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_d2backward = (float*)malloc(sizeof(float) * 4 * p *p);

    float* leg_trans      = (float*)malloc(sizeof(float) * leg_mat_size);
    float* leg_trans_inv  = (float*)malloc(sizeof(float) * leg_mat_size);
    float* d1_leg_trans   = (float*)malloc(sizeof(float) * leg_mat_size);
    float* d2_leg_trans   = (float*)malloc(sizeof(float) * leg_mat_size);

    BlasSht c(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt", dft_forward,
        dft_backward, dft_d1backward, dft_d2backward, leg_trans,
        leg_trans_inv, d1_leg_trans, d2_leg_trans);
    //c.test(num_vesicles);

    int vec_length = 2 * p * (p + 1) * num_vesicles;
    scalar *inputs, *outputs, *outputs_2, *work_arr;
    inputs    = (scalar*) malloc(vec_length * sizeof(scalar));
    outputs   = (scalar*) malloc(vec_length * sizeof(scalar));
    outputs_2 = (scalar*) malloc(vec_length * sizeof(scalar));
    work_arr  = (scalar*) malloc(2 * vec_length * sizeof(scalar));
    for (int i = 0; i < vec_length; i++) {
        inputs[i] = (float) rand() / (float) RAND_MAX;
    }
    c.forward(inputs, work_arr, num_vesicles, outputs);
    c.backward_du(outputs, work_arr, num_vesicles, outputs_2);

    //Testing the device
    DeviceCPU<float> cpu1;
    cpu1.InitializeSHT(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt");
  
    scalar *du, *dv, *shc;
    du  = (scalar*) malloc(vec_length * sizeof(scalar));
    dv  = (scalar*) malloc(vec_length * sizeof(scalar));
    shc = (scalar*) malloc(vec_length * sizeof(scalar));
    
    cpu1.ShAna(inputs, work_arr, num_vesicles, outputs);
    cpu1.ShSynDu(outputs, work_arr, num_vesicles, du);
    for (int i = 0; i < 2 * p * (p + 1) * num_vesicles; i++) {
        cout<<outputs_2[i]<<" "<<du[i]-outputs_2[i]<<endl;
    }
    
    Scalars<float> sc(cpu1, p, num_vesicles,inputs);
    Scalars<float> dsc(cpu1, p, num_vesicles);
    Scalars<float> trash(cpu1, p, num_vesicles);
    SHTrans<float> diff(cpu1, p, num_vesicles);
    
    diff.FirstDerivatives(sc,dsc,trash);
    
    for (int i = 0; i < 2 * p * (p + 1) * num_vesicles; i++) {
        cout<<outputs_2[i]<<" "<<dsc.data_[i]-du[i]<<endl;
    }
    
    
    free(dft_forward);
    free(dft_backward);
    free(dft_d1backward);
    free(dft_d2backward);
    free(leg_trans);
    free(leg_trans_inv);
    free(d1_leg_trans);
    free(d2_leg_trans);
    free(inputs);
    free(outputs);
    free(outputs_2);
    free(work_arr);
    free(du);
    free(dv);
    free(shc);

    return 0;
}
