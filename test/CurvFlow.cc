#include<exception>
#include "SHVectors.h"
#include "Surface.h"
#include<iostream>

using namespace std;

extern void CurvFlow(float* x,int vec_len,float ts,int m,float *x_final)
{
    int p = 12;
    int np = vec_len/3;
    SHVectors<float> pos_vec(p,1,x);

    // initializing Vesicle
    Surface<float> S(p,1,pos_vec);

    SHVectors<float> bend_force(p,1);

    for(int tt=0;tt<m;++tt)
    {
        AxPy(S.k_, S.normal_, (float) 0.0, bend_force);
        //Add pressure to preserve volume
        AxPy(-ts,bend_force, S.x_, pos_vec);
        S.SetX(pos_vec);
        cout<<tt<<" ";
    }
    cout<<endl;
    
    for(int idx=0;idx<np;++idx)
        x_final[idx] = S.x_.data_[idx];
    
    return;
}
