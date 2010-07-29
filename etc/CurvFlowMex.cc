#include "mex.h"

void CurvFlow(float* x,int vec_len,float ts,int m,float *x_final);

//Matlab syntax Xf = CurvFlow(X, ts, m);
void mexFunction(
    int             nlhs,
    mxArray      *plhs[],
    int             nrhs,
    const mxArray *prhs[])
{
    if (nrhs != 3) {
        mexErrMsgTxt("CurvFlowMex requires three input arguments.");
    } else if (nlhs > 1) {
        mexErrMsgTxt("CurvFlowMex requires one output argument.");
    }

    int np, m;
    double *XX_d, *ts_d, *m_d;
    float ts, *x, *x_final;

    XX_d = (double *) mxGetPr(prhs[0]);
    ts_d = (double *) mxGetPr(prhs[1]);
    m_d  = (double *) mxGetPr(prhs[2]);
    
    m = (int) *m_d;
    ts = (float) *ts_d;

    np = mxGetM(prhs[0]);
    x = new float[np];
    x_final = new float[np];

    for(int idx=0;idx<np;++idx)
        x[idx] = (float) XX_d[idx];

    CurvFlow(x,np,ts,m,x_final);
    
    //output
    plhs[0] = mxCreateDoubleMatrix(np,1,mxREAL);
    double *x_out = mxGetPr(plhs[0]);

    for(int idx=0;idx<np;++idx)
        x_out[idx] = (double) x_final[idx];

    delete []x;
    delete []x_final;
    
    return;
}
