#include<iostream>
#include "DeviceCPU.h"
#include "Scalars.h"
#include "Vectors.h"
#include "SHTrans.h"


int main(int argc, char ** argv)
{
    const int p = 12;
    const int num_vesicles = 3;
    int leg_mat_size      = (p + 1) * (p + 1) * (p + 2);
    float* dft_forward    = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_backward   = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_d1backward = (float*)malloc(sizeof(float) * 4 * p *p);
    float* dft_d2backward = (float*)malloc(sizeof(float) * 4 * p *p);

    float* leg_trans      = (float*)malloc(sizeof(float) * leg_mat_size);
    float* leg_trans_inv  = (float*)malloc(sizeof(float) * leg_mat_size);
    float* d1_leg_trans   = (float*)malloc(sizeof(float) * leg_mat_size);
    float* d2_leg_trans   = (float*)malloc(sizeof(float) * leg_mat_size);

    float* trans_in       = (float*)malloc(sizeof(float) * 2 * p * (p + 1) * num_vesicles);
    float* trans_out      = (float*)malloc(sizeof(float) * 2 * p * (p + 1) * num_vesicles);

    BlasSht c(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt", dft_forward,
        dft_backward, dft_d1backward, dft_d2backward, leg_trans,
        leg_trans_inv, d1_leg_trans, d2_leg_trans, trans_in,
        trans_out);
    c.test(num_vesicles);

    //
    int vec_length = 2 * p * (p + 1) * num_vesicles;
    scalar *inputs, *outputs, *outputs_2;
    inputs    = (scalar*) malloc(vec_length * sizeof(scalar));
    outputs   = (scalar*) malloc(vec_length * sizeof(scalar));
    outputs_2 = (scalar*) malloc(vec_length * sizeof(scalar));

    for (int i = 0; i < vec_length; i++) {
        inputs[i] = (float) rand() / (float) RAND_MAX;
    }
    c.forward(inputs, num_vesicles, outputs);
    c.backward_du(outputs, num_vesicles, outputs_2);

    DeviceCPU<float> cpu1;
//     cpu1.sht_.InitializeBlasSht(p, "../data/legTrans12_single.txt",
//         "../data/legTransInv12_single.txt",
//         "../data/d1legTrans12_single.txt",
//         "../data/d2legTrans12_single.txt", dft_forward,
//         dft_backward, dft_d1backward, dft_d2backward, leg_trans,
//         leg_trans_inv, d1_leg_trans, d2_leg_trans, trans_in,
//         trans_out);
    
    cpu1.trans_in       = trans_in;      
    cpu1.trans_out      = trans_out;

    cpu1.InitializeSHT(p, "../data/legTrans12_single.txt",
        "../data/legTransInv12_single.txt",
        "../data/d1legTrans12_single.txt",
        "../data/d2legTrans12_single.txt");
  
    scalar *du, *dv, *shc;
    du  = (scalar*) malloc(vec_length * sizeof(scalar));
    dv  = (scalar*) malloc(vec_length * sizeof(scalar));
    shc = (scalar*) malloc(vec_length * sizeof(scalar));
    
    cpu1.ShAna(inputs, num_vesicles, outputs);
    cpu1.ShSynDu(outputs, num_vesicles, du);
    for (int i = 0; i < 2 * p * (p + 1) * num_vesicles; i++) {
        cout<<outputs_2[i]<<" "<<du[i]-outputs_2[i]<<endl;
    }
    
    //   Scalars<float> sc(cpu1, p, num_vesicles,inputs);
    //   Scalars<float> dsc(cpu1, p, num_vesicles);
    //   Scalars<float> trash(cpu1, p, num_vesicles);
    
    //   SHTrans<float> diff(cpu1, p, num_vesicles);
    
    //   diff.AllDerivatives(sc,trash,trash,trash,trash,dsc);
    //   //diff.FirstDerivatives(sc,dusc,dvsc);
    
    //   for (int i = 0; i < 2 * p * (p + 1) * num_vesicles; i++) {
    //       cout<<outputs_2[i]<<" "<<dsc.data_[i]-outputs_2[i]<<endl;
    //   }
    
    
    free(dft_forward);
    free(dft_backward);
    free(dft_d1backward);
    free(dft_d2backward);
    free(leg_trans);
    free(leg_trans_inv);
    free(d1_leg_trans);
    free(d2_leg_trans);
    free(trans_in);
    free(trans_out);
    free(inputs);
    free(outputs);
    free(outputs_2);
    free(du);
    free(dv);
    free(shc);

    return 0;
}
