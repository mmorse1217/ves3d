#include "CPUKernels.h"

#define I_PI 1.0/M_PI/8.0
#define IDEAL_ALIGNMENT 16
#define SIMD_LEN_F (IDEAL_ALIGNMENT / sizeof(float))
#define SIMD_LEN_D (IDEAL_ALIGNMENT / sizeof(double))

void DirectStokesKernel(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *qw, const float *trg, const float *src,
    const float *den, float *pot)
{
    float tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    
#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for (int vt=0; vt<n_surfs; vt++)
    {
        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            px = 0;
            py = 0;
            pz = 0;
            
            tx=trg[3*vt*stride +                   trg_idx];
            ty=trg[3*vt*stride + stride +          trg_idx];
            tz=trg[3*vt*stride + stride + stride + trg_idx];
            
            for (int s=0; s<stride; s++)
            {
                dx=src[3*stride*vt +                   s]-tx;
                dy=src[3*stride*vt + stride +          s]-ty;
                dz=src[3*stride*vt + stride + stride + s]-tz;

                invR = dx*dx;
                invR+= dy*dy;
                invR+= dz*dz;
                
                if (invR!=0)
                    invR = 1.0/sqrt(invR);
            
                cpx = den[3*stride*vt +                   s] * qw[s]; 
                cpy = den[3*stride*vt + stride +          s] * qw[s]; 
                cpz = den[3*stride*vt + stride + stride + s] * qw[s]; 
                
                cc  = dx*cpx;
                cc += dy*cpy;
                cc += dz*cpz;
                cc *= invR;
                cc *= invR;

                cpx += cc*dx;
                cpy += cc*dy;
                cpz += cc*dz;
                
                px += cpx*invR;
                py += cpy*invR;
                pz += cpz*invR;
            }
            pot[3*vt*stride +                  trg_idx] = px * I_PI;
            pot[3*vt*stride + stride +         trg_idx] = py * I_PI;
            pot[3*vt*stride + stride +stride + trg_idx] = pz * I_PI;
        }

    }
}

void DirectStokesKernel_Noqw(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *trg, const float *src, 
    const float *den, float *pot)
{
    float tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;

#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for (int vt=0; vt<n_surfs; vt++)
    {
        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            px = 0;
            py = 0;
            pz = 0;
            
            tx=trg[3*vt*stride +                   trg_idx];
            ty=trg[3*vt*stride + stride +          trg_idx];
            tz=trg[3*vt*stride + stride + stride + trg_idx];
            
            for (int s=0; s<stride; s++)
            {
                dx=src[3*stride*vt +                   s]-tx;
                dy=src[3*stride*vt + stride +          s]-ty;
                dz=src[3*stride*vt + stride + stride + s]-tz;

                invR = dx*dx;
                invR+= dy*dy;
                invR+= dz*dz;
                
                if (invR!=0)
                    invR = 1.0/sqrt(invR);
            
                cpx = den[3*stride*vt +                   s]; 
                cpy = den[3*stride*vt + stride +          s]; 
                cpz = den[3*stride*vt + stride + stride + s]; 
                
                cc  = dx*cpx;
                cc += dy*cpy;
                cc += dz*cpz;
                cc *= invR;
                cc *= invR;

                cpx += cc*dx;
                cpy += cc*dy;
                cpz += cc*dz;
                
                px += cpx*invR;
                py += cpy*invR;
                pz += cpz*invR;
            }
            pot[3*vt*stride +                  trg_idx] = px * I_PI;
            pot[3*vt*stride + stride +         trg_idx] = py * I_PI;
            pot[3*vt*stride + stride +stride + trg_idx] = pz * I_PI;
        }
    }
}

void DirectStokesKernel(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *qw, const double *trg, 
    const double *src, const double *den, double *pot)
{
    double tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    
#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for (int vt=0; vt<n_surfs; vt++)
    {
        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            px = 0;
            py = 0;
            pz = 0;
            
            tx=trg[3*vt*stride +                   trg_idx];
            ty=trg[3*vt*stride + stride +          trg_idx];
            tz=trg[3*vt*stride + stride + stride + trg_idx];
            
            for (int s=0; s<stride; s++)
            {
                dx=src[3*stride*vt +                   s]-tx;
                dy=src[3*stride*vt + stride +          s]-ty;
                dz=src[3*stride*vt + stride + stride + s]-tz;

                invR = dx*dx;
                invR+= dy*dy;
                invR+= dz*dz;
                
                if (invR!=0)
                    invR = 1.0/sqrt(invR);
            
                cpx = den[3*stride*vt +                   s] * qw[s]; 
                cpy = den[3*stride*vt + stride +          s] * qw[s]; 
                cpz = den[3*stride*vt + stride + stride + s] * qw[s]; 
                
                cc  = dx*cpx;
                cc += dy*cpy;
                cc += dz*cpz;
                cc *= invR;
                cc *= invR;

                cpx += cc*dx;
                cpy += cc*dy;
                cpz += cc*dz;
                
                px += cpx*invR;
                py += cpy*invR;
                pz += cpz*invR;
            }
            pot[3*vt*stride +                  trg_idx] = px * I_PI;
            pot[3*vt*stride + stride +         trg_idx] = py * I_PI;
            pot[3*vt*stride + stride +stride + trg_idx] = pz * I_PI;
        }
    }
}

void DirectStokesKernel_Noqw(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const double *trg, const double *src, 
    const double *den, double *pot)
{

    double tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;

#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
    for (int vt=0; vt<n_surfs; vt++)
    {
        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            px = 0;
            py = 0;
            pz = 0;
            
            tx=trg[3*vt*stride +                   trg_idx];
            ty=trg[3*vt*stride + stride +          trg_idx];
            tz=trg[3*vt*stride + stride + stride + trg_idx];
            
            for (int s=0; s<stride; s++)
            {
                dx=src[3*stride*vt +                   s]-tx;
                dy=src[3*stride*vt + stride +          s]-ty;
                dz=src[3*stride*vt + stride + stride + s]-tz;

                invR = dx*dx;
                invR+= dy*dy;
                invR+= dz*dz;
                
                if (invR!=0)
                    invR = 1.0/sqrt(invR);
            
                cpx = den[3*stride*vt +                   s]; 
                cpy = den[3*stride*vt + stride +          s]; 
                cpz = den[3*stride*vt + stride + stride + s]; 
                
                cc  = dx*cpx;
                cc += dy*cpy;
                cc += dz*cpz;
                cc *= invR;
                cc *= invR;

                cpx += cc*dx;
                cpy += cc*dy;
                cpz += cc*dz;
                
                px += cpx*invR;
                py += cpy*invR;
                pz += cpz*invR;
            }
            pot[3*vt*stride +                  trg_idx] = px * I_PI;
            pot[3*vt*stride + stride +         trg_idx] = py * I_PI;
            pot[3*vt*stride + stride +stride + trg_idx] = pz * I_PI;
        }
    }
} 

///////////////////////////////////////////////////////////////////////////////////
void DirectStokesSSE(int stride, int n_surfs, int trg_idx_head, int trg_idx_tail, 
    const float *qw, const float *trg, const float *src, 
    const float *den, float *pot)
{
    if (stride%4) // necessary for proper alignment of sources
        exit(1);

#pragma omp parallel for
    ///@todo add the openmp instructions
    for (int vt=0; vt<n_surfs; vt++)
    {
        float aux_arr[3*SIMD_LEN_F+3]; 
        float *tempvalx; 
        float *tempvaly; 
        float *tempvalz; 
        size_t residual = size_t(aux_arr)%IDEAL_ALIGNMENT;
        if (residual)  // if aux_arr is misaligned
            tempvalx = aux_arr + (IDEAL_ALIGNMENT - residual)/sizeof(float);
        else tempvalx = aux_arr;
        if (size_t(tempvalx)%IDEAL_ALIGNMENT)  // for debugging
            exit(1);
        tempvaly=tempvalx+SIMD_LEN_F;
        tempvalz=tempvaly+SIMD_LEN_F;

        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            float p[3]={0,0,0};
            float tx=trg[3*vt*stride            + trg_idx];
            float ty=trg[3*vt*stride +   stride + trg_idx];
            float tz=trg[3*vt*stride + 2*stride + trg_idx];

            residual = size_t(src+ 3*stride*vt)%IDEAL_ALIGNMENT;
            if (residual)
                residual = IDEAL_ALIGNMENT - residual;

            // Handle start data if it is not 16-byte aligned
            size_t s;
            // residual = stride;
            for (s=0; s<residual; s++)
            {
                float dX_reg=src[3*stride*vt+           s]-tx;
                float dY_reg=src[3*stride*vt+  stride + s]-ty;
                float dZ_reg=src[3*stride*vt+2*stride + s]-tz;

                float invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
                if (invR!=0)
                    invR = 1.0/sqrt(invR);

                float cur_pot_x = den[3*stride*vt +           s] * qw[s];
                float cur_pot_y = den[3*stride*vt +  stride + s] * qw[s];
                float cur_pot_z = den[3*stride*vt +2*stride + s] * qw[s];

                float tmp_scalar = (dX_reg*cur_pot_x + dY_reg*cur_pot_y + 
                    dZ_reg*cur_pot_z)*invR*invR;
                
                cur_pot_x += tmp_scalar*dX_reg;
                cur_pot_y += tmp_scalar*dY_reg;
                cur_pot_z += tmp_scalar*dZ_reg;

                p[0] += cur_pot_x*invR;
                p[1] += cur_pot_y*invR;
                p[2] += cur_pot_z*invR;
            }

            __m128 txi = _mm_load1_ps (&tx);
            __m128 tyi = _mm_load1_ps (&ty);
            __m128 tzi = _mm_load1_ps (&tz);

            // for (int s=0; s<stride; s++)
            __m128 tempx;
            __m128 tempy;
            __m128 tempz;

            tempx = _mm_setzero_ps();
            tempy = _mm_setzero_ps();
            tempz = _mm_setzero_ps();

            // Load and calculate in groups of SIMD_LEN
            size_t loop_limit = stride-SIMD_LEN_F;
            for (; s <= loop_limit; s += SIMD_LEN_F) {
                __m128 sxj = _mm_load_ps (src+3*stride*vt+         s);
                __m128 syj = _mm_load_ps (src+3*stride*vt+  stride+s);
                __m128 szj = _mm_load_ps (src+3*stride*vt+2*stride+s);

                // this could be vectorized assuming den and q are 16-byte aligned
                __m128 sdenx = _mm_set_ps (
                    den[3*stride*vt +s+3] * qw[s+3],
                    den[3*stride*vt +s+2] * qw[s+2],
                    den[3*stride*vt +s+1] * qw[s+1],
                    den[3*stride*vt +s] * qw[s]);

                __m128 sdeny = _mm_set_ps (
                    den[3*stride*vt+stride +s+3] * qw[s+3],
                    den[3*stride*vt+stride +s+2] * qw[s+2],
                    den[3*stride*vt+stride +s+1] * qw[s+1],
                    den[3*stride*vt+stride +s] * qw[s]);

                __m128 sdenz = _mm_set_ps (
                    den[3*stride*vt+2*stride +s+3] * qw[s+3],
                    den[3*stride*vt+2*stride +s+2] * qw[s+2],
                    den[3*stride*vt+2*stride +s+1] * qw[s+1],
                    den[3*stride*vt+2*stride +s] * qw[s]
                                           );

                //__m128 sdenx = _mm_load_ps (src+3*stride*vt+         s);
                //__m128 sdeny = _mm_load_ps (src+3*stride*vt+  stride+s);
                //__m128 sdenz = _mm_load_ps (src+3*stride*vt+2*stride+s);

                __m128 dX, dY, dZ;
                __m128 dR2;
                __m128 S;

                dX = _mm_sub_ps(txi , sxj);
                dY = _mm_sub_ps(tyi , syj);
                dZ = _mm_sub_ps(tzi , szj);

                sxj = _mm_mul_ps(dX, dX); 
                syj = _mm_mul_ps(dY, dY);
                szj = _mm_mul_ps(dZ, dZ);

                dR2 = _mm_add_ps(sxj, syj);
                dR2 = _mm_add_ps(szj, dR2);

                __m128 zero = _mm_setzero_ps ();
                __m128 is_zero = _mm_cmpeq_ps(dR2, zero);

                // S = _mm_rsqrt_ps(dR2);
                const __m128 approx = _mm_rsqrt_ps( dR2 );
                const __m128 muls = _mm_mul_ps(_mm_mul_ps(dR2, approx), approx);
                const __m128 three = _mm_set1_ps (3.0f);
                const __m128 half4 = _mm_set1_ps (0.5f);
                S = _mm_mul_ps(_mm_mul_ps(half4, approx), _mm_sub_ps(three, muls));
                S = _mm_andnot_ps (is_zero, S);

                __m128 dotx = _mm_mul_ps (dX, sdenx);
                __m128 doty = _mm_mul_ps (dY, sdeny);
                __m128 dotz = _mm_mul_ps (dZ, sdenz);

                __m128 dot_sum = _mm_add_ps (dotx, doty);
                dot_sum = _mm_add_ps (dot_sum, dotz);

                dot_sum = _mm_mul_ps (dot_sum, S);
                dot_sum = _mm_mul_ps (dot_sum, S);

                dotx = _mm_mul_ps (dot_sum, dX);
                doty = _mm_mul_ps (dot_sum, dY);
                dotz = _mm_mul_ps (dot_sum, dZ);

                sdenx = _mm_add_ps (sdenx, dotx);
                sdeny = _mm_add_ps (sdeny, doty);
                sdenz = _mm_add_ps (sdenz, dotz);

                sdenx = _mm_mul_ps (sdenx, S);
                sdeny = _mm_mul_ps (sdeny, S);
                sdenz = _mm_mul_ps (sdenz, S);

                tempx = _mm_add_ps (sdenx, tempx);
                tempy = _mm_add_ps (sdeny, tempy);
                tempz = _mm_add_ps (sdenz, tempz);

            }

            _mm_store_ps(tempvalx, tempx); 
            _mm_store_ps(tempvaly, tempy); 
            _mm_store_ps(tempvalz, tempz); 

            for (size_t k = 0; k < SIMD_LEN_F; k++) {
                p[0] += tempvalx[k];
                p[1] += tempvaly[k];
                p[2] += tempvalz[k];
            }

            if (s!=size_t(stride))
                exit(1);

            pot[3*vt*stride +            trg_idx] = p[0] * I_PI;
            pot[3*vt*stride +   stride + trg_idx] = p[1] * I_PI;
            pot[3*vt*stride + 2*stride + trg_idx] = p[2] * I_PI;
        }
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////////
void DirectStokesSSE_Noqw(int stride, int n_surfs, int trg_idx_head,
    int trg_idx_tail, const float *trg, const float *src, 
    const float *den, float *pot)
{
    if (stride%4) // necessary for proper alignment of sources
        exit(1);
    
#pragma omp parallel for
    ///@todo add the openmp instructions
    for (int vt=0; vt<n_surfs; vt++)
    {
        float aux_arr[3*SIMD_LEN_F+3]; 
        float *tempvalx; 
        float *tempvaly; 
        float *tempvalz; 
        size_t residual = size_t(aux_arr)%IDEAL_ALIGNMENT;
        if (residual)  // if aux_arr is misaligned
            tempvalx = aux_arr + (IDEAL_ALIGNMENT - residual)/sizeof(float);
        else tempvalx = aux_arr;
        if (size_t(tempvalx)%IDEAL_ALIGNMENT)  // for debugging
            exit(1);
        tempvaly=tempvalx+SIMD_LEN_F;
        tempvalz=tempvaly+SIMD_LEN_F;

        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            float p[3]={0,0,0};
            float tx=trg[3*vt*stride            + trg_idx];
            float ty=trg[3*vt*stride +   stride + trg_idx];
            float tz=trg[3*vt*stride + 2*stride + trg_idx];

            residual = size_t(src+ 3*stride*vt)%IDEAL_ALIGNMENT;
            if (residual)
                residual = IDEAL_ALIGNMENT - residual;

            // Handle start data if it is not 16-byte aligned
            size_t s;
            // residual = stride;
            for (s=0; s<residual; s++)
            {
                float dX_reg=src[3*stride*vt+           s]-tx;
                float dY_reg=src[3*stride*vt+  stride + s]-ty;
                float dZ_reg=src[3*stride*vt+2*stride + s]-tz;

                float invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
                if (invR!=0)
                    invR = 1.0/sqrt(invR);

                float cur_pot_x = den[3*stride*vt +           s];
                float cur_pot_y = den[3*stride*vt +  stride + s];
                float cur_pot_z = den[3*stride*vt +2*stride + s];

                float tmp_scalar = (dX_reg*cur_pot_x + dY_reg*cur_pot_y 
                    + dZ_reg*cur_pot_z)*invR*invR;
                cur_pot_x += tmp_scalar*dX_reg;
                cur_pot_y += tmp_scalar*dY_reg;
                cur_pot_z += tmp_scalar*dZ_reg;

                p[0] += cur_pot_x*invR;
                p[1] += cur_pot_y*invR;
                p[2] += cur_pot_z*invR;
            }

            __m128 txi = _mm_load1_ps (&tx);
            __m128 tyi = _mm_load1_ps (&ty);
            __m128 tzi = _mm_load1_ps (&tz);

            // for (int s=0; s<stride; s++)
            __m128 tempx;
            __m128 tempy;
            __m128 tempz;

            tempx = _mm_setzero_ps();
            tempy = _mm_setzero_ps();
            tempz = _mm_setzero_ps();

            // Load and calculate in groups of SIMD_LEN
            size_t loop_limit = stride-SIMD_LEN_F;
            for (; s <= loop_limit; s += SIMD_LEN_F) {
                __m128 sxj = _mm_load_ps (src+3*stride*vt+         s);
                __m128 syj = _mm_load_ps (src+3*stride*vt+  stride+s);
                __m128 szj = _mm_load_ps (src+3*stride*vt+2*stride+s);

                // this could be vectorized assuming den and q are 16-byte aligned
                __m128 sdenx = _mm_set_ps (
                    den[3*stride*vt +s+3],
                    den[3*stride*vt +s+2],
                    den[3*stride*vt +s+1],
                    den[3*stride*vt +s]);

                __m128 sdeny = _mm_set_ps (
                    den[3*stride*vt+stride +s+3],
                    den[3*stride*vt+stride +s+2],
                    den[3*stride*vt+stride +s+1],
                    den[3*stride*vt+stride +s]);

                __m128 sdenz = _mm_set_ps (
                    den[3*stride*vt+2*stride +s+3],
                    den[3*stride*vt+2*stride +s+2],
                    den[3*stride*vt+2*stride +s+1],
                    den[3*stride*vt+2*stride +s]);

                //__m128 sdenx = _mm_load_ps (src+3*stride*vt+         s);
                //__m128 sdeny = _mm_load_ps (src+3*stride*vt+  stride+s);
                //__m128 sdenz = _mm_load_ps (src+3*stride*vt+2*stride+s);

                __m128 dX, dY, dZ;
                __m128 dR2;
                __m128 S;

                dX = _mm_sub_ps(txi , sxj);
                dY = _mm_sub_ps(tyi , syj);
                dZ = _mm_sub_ps(tzi , szj);

                sxj = _mm_mul_ps(dX, dX); 
                syj = _mm_mul_ps(dY, dY);
                szj = _mm_mul_ps(dZ, dZ);

                dR2 = _mm_add_ps(sxj, syj);
                dR2 = _mm_add_ps(szj, dR2);

                __m128 zero = _mm_setzero_ps ();
                __m128 is_zero = _mm_cmpeq_ps(dR2, zero);

                // S = _mm_rsqrt_ps(dR2);
                const __m128 approx = _mm_rsqrt_ps( dR2 );
                const __m128 muls = _mm_mul_ps(_mm_mul_ps(dR2, approx), approx);
                const __m128 three = _mm_set1_ps (3.0f);
                const __m128 half4 = _mm_set1_ps (0.5f);
                S = _mm_mul_ps(_mm_mul_ps(half4, approx), _mm_sub_ps(three, muls));
                S = _mm_andnot_ps (is_zero, S);

                __m128 dotx = _mm_mul_ps (dX, sdenx);
                __m128 doty = _mm_mul_ps (dY, sdeny);
                __m128 dotz = _mm_mul_ps (dZ, sdenz);

                __m128 dot_sum = _mm_add_ps (dotx, doty);
                dot_sum = _mm_add_ps (dot_sum, dotz);

                dot_sum = _mm_mul_ps (dot_sum, S);
                dot_sum = _mm_mul_ps (dot_sum, S);

                dotx = _mm_mul_ps (dot_sum, dX);
                doty = _mm_mul_ps (dot_sum, dY);
                dotz = _mm_mul_ps (dot_sum, dZ);

                sdenx = _mm_add_ps (sdenx, dotx);
                sdeny = _mm_add_ps (sdeny, doty);
                sdenz = _mm_add_ps (sdenz, dotz);

                sdenx = _mm_mul_ps (sdenx, S);
                sdeny = _mm_mul_ps (sdeny, S);
                sdenz = _mm_mul_ps (sdenz, S);

                tempx = _mm_add_ps (sdenx, tempx);
                tempy = _mm_add_ps (sdeny, tempy);
                tempz = _mm_add_ps (sdenz, tempz);

            }

            _mm_store_ps(tempvalx, tempx); 
            _mm_store_ps(tempvaly, tempy); 
            _mm_store_ps(tempvalz, tempz); 

            for (size_t k = 0; k < SIMD_LEN_F; k++) {
                p[0] += tempvalx[k];
                p[1] += tempvaly[k];
                p[2] += tempvalz[k];
            }

            if (s!=size_t(stride))
                exit(1);

            pot[3*vt*stride +            trg_idx] = p[0] * I_PI;
            pot[3*vt*stride +   stride + trg_idx] = p[1] * I_PI;
            pot[3*vt*stride + 2*stride + trg_idx] = p[2] * I_PI;
        }
    }

    return;
}

////////////////////////////////////////////////////////////////////////////
void DirectStokesSSE(int stride, int n_surfs, int trg_idx_head, 
    int trg_idx_tail, const double *qw, const double *trg, 
    const double *src, const double *den, double *pot)
{
    if (stride%4) // necessary for proper alignment of sources
        exit(1);
    
#pragma omp parallel for
    ///@todo add the openmp instructions
    for (int vt=0; vt<n_surfs; vt++)
    {
        double aux_arr[3*SIMD_LEN_D+3]; 
        double *tempvalx; 
        double *tempvaly; 
        double *tempvalz; 
        size_t residual = size_t(aux_arr)%IDEAL_ALIGNMENT;
        if (residual)  // if aux_arr is misaligned
            tempvalx = aux_arr + (IDEAL_ALIGNMENT - residual)/sizeof(double);
        else tempvalx = aux_arr;
        if (size_t(tempvalx)%IDEAL_ALIGNMENT)  // for debugging
            exit(1);
        tempvaly=tempvalx+SIMD_LEN_D;
        tempvalz=tempvaly+SIMD_LEN_D;

        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            double p[3]={0,0,0};
            double tx=trg[3*vt*stride            + trg_idx];
            double ty=trg[3*vt*stride +   stride + trg_idx];
            double tz=trg[3*vt*stride + 2*stride + trg_idx];

            residual = size_t(src+ 3*stride*vt)%IDEAL_ALIGNMENT;
            if (residual)
                residual = IDEAL_ALIGNMENT - residual;

            // Handle start data if it is not 16-byte aligned
            size_t s;
            // residual = stride;
            for (s=0; s<residual; s++)
            {
                double dX_reg=src[3*stride*vt+           s]-tx;
                double dY_reg=src[3*stride*vt+  stride + s]-ty;
                double dZ_reg=src[3*stride*vt+2*stride + s]-tz;

                double invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
                if (invR!=0)
                    invR = 1.0/sqrt(invR);

                double cur_pot_x = den[3*stride*vt +           s] * qw[s];
                double cur_pot_y = den[3*stride*vt +  stride + s] * qw[s];
                double cur_pot_z = den[3*stride*vt +2*stride + s] * qw[s];

                double tmp_scalar = (dX_reg*cur_pot_x + dY_reg*cur_pot_y
                    + dZ_reg*cur_pot_z)*invR*invR;
                
                cur_pot_x += tmp_scalar*dX_reg;
                cur_pot_y += tmp_scalar*dY_reg;
                cur_pot_z += tmp_scalar*dZ_reg;

                p[0] += cur_pot_x*invR;
                p[1] += cur_pot_y*invR;
                p[2] += cur_pot_z*invR;
            }

            __m128d txi = _mm_load1_pd (&tx);
            __m128d tyi = _mm_load1_pd (&ty);
            __m128d tzi = _mm_load1_pd (&tz);

            // for (int s=0; s<stride; s++)
            __m128d tempx;
            __m128d tempy;
            __m128d tempz;

            tempx = _mm_setzero_pd();
            tempy = _mm_setzero_pd();
            tempz = _mm_setzero_pd();

            // Load and calculate in groups of SIMD_LEN
            size_t loop_limit = stride-SIMD_LEN_D;
            for (; s <= loop_limit; s += SIMD_LEN_D) {
                __m128d sxj = _mm_load_pd (src+3*stride*vt+         s);
                __m128d syj = _mm_load_pd (src+3*stride*vt+  stride+s);
                __m128d szj = _mm_load_pd (src+3*stride*vt+2*stride+s);

                // this could be vectorized assuming den and q are 16-byte aligned
                __m128d sdenx = _mm_set_pd(
                    den[3*stride*vt +s+1] * qw[s+1],
                    den[3*stride*vt +s] * qw[s]);

                __m128d sdeny = _mm_set_pd(
                    den[3*stride*vt+stride +s+1] * qw[s+1],
                    den[3*stride*vt+stride +s] * qw[s]);

                __m128d sdenz = _mm_set_pd (
                    den[3*stride*vt+2*stride +s+1] * qw[s+1],
                    den[3*stride*vt+2*stride +s] * qw[s]);

                __m128d dX, dY, dZ;
                __m128d dR2;
                __m128d S;

                dX = _mm_sub_pd(txi , sxj);
                dY = _mm_sub_pd(tyi , syj);
                dZ = _mm_sub_pd(tzi , szj);

                sxj = _mm_mul_pd(dX, dX); 
                syj = _mm_mul_pd(dY, dY);
                szj = _mm_mul_pd(dZ, dZ);

                dR2 = _mm_add_pd(sxj, syj);
                dR2 = _mm_add_pd(szj, dR2);

                __m128d zero = _mm_setzero_pd ();
                __m128d is_zero = _mm_cmpeq_pd(dR2, zero);
                
                ///@bug this part of code may be buggy -abtin
                // S = _mm_rsqrt_ps(dR2);
                __m128 dR2_s  = _mm_cvtpd_ps(dR2);
                
                const __m128 approx = _mm_rsqrt_ps( dR2_s );
                const __m128 muls  = _mm_mul_ps(_mm_mul_ps(dR2_s, approx), approx);
                const __m128 three = _mm_set1_ps (3.0f);
                const __m128 half4 = _mm_set1_ps (0.5f);
                __m128 S_s =  _mm_mul_ps(_mm_mul_ps(half4, approx), _mm_sub_ps(three, muls));
                S = _mm_cvtps_pd(S_s);
                S = _mm_andnot_pd (is_zero, S);
                ///@bug

                __m128d dotx = _mm_mul_pd (dX, sdenx);
                __m128d doty = _mm_mul_pd (dY, sdeny);
                __m128d dotz = _mm_mul_pd (dZ, sdenz);

                __m128d dot_sum = _mm_add_pd (dotx, doty);
                dot_sum = _mm_add_pd (dot_sum, dotz);

                dot_sum = _mm_mul_pd (dot_sum, S);
                dot_sum = _mm_mul_pd (dot_sum, S);

                dotx = _mm_mul_pd (dot_sum, dX);
                doty = _mm_mul_pd (dot_sum, dY);
                dotz = _mm_mul_pd (dot_sum, dZ);

                sdenx = _mm_add_pd (sdenx, dotx);
                sdeny = _mm_add_pd (sdeny, doty);
                sdenz = _mm_add_pd (sdenz, dotz);

                sdenx = _mm_mul_pd (sdenx, S);
                sdeny = _mm_mul_pd (sdeny, S);
                sdenz = _mm_mul_pd (sdenz, S);

                tempx = _mm_add_pd (sdenx, tempx);
                tempy = _mm_add_pd (sdeny, tempy);
                tempz = _mm_add_pd (sdenz, tempz);

            }

            _mm_store_pd(tempvalx, tempx); 
            _mm_store_pd(tempvaly, tempy); 
            _mm_store_pd(tempvalz, tempz); 

            for (size_t k = 0; k < SIMD_LEN_D; k++) {
                p[0] += tempvalx[k];
                p[1] += tempvaly[k];
                p[2] += tempvalz[k];
            }

            if (s!=size_t(stride))
                exit(1);
            
            pot[3*vt*stride +            trg_idx] = p[0] * I_PI;
            pot[3*vt*stride +   stride + trg_idx] = p[1] * I_PI;
            pot[3*vt*stride + 2*stride + trg_idx] = p[2] * I_PI;
        }
    }

    return;
}

void DirectStokesSSE_Noqw(int stride, int n_surfs, int trg_idx_head, 
    int trg_idx_tail, const double *trg, const double *src,
    const double *den, double *pot)
{
    if (stride%4) // necessary for proper alignment of sources
        exit(1);
    
#pragma omp parallel for
    ///@todo add the openmp instructions
    for (int vt=0; vt<n_surfs; vt++)
    {
        double aux_arr[3*SIMD_LEN_D+3]; 
        double *tempvalx; 
        double *tempvaly; 
        double *tempvalz; 
        size_t residual = size_t(aux_arr)%IDEAL_ALIGNMENT;
        if (residual)  // if aux_arr is misaligned
            tempvalx = aux_arr + (IDEAL_ALIGNMENT - residual)/sizeof(double);
        else tempvalx = aux_arr;
        if (size_t(tempvalx)%IDEAL_ALIGNMENT)  // for debugging
            exit(1);
        tempvaly=tempvalx+SIMD_LEN_D;
        tempvalz=tempvaly+SIMD_LEN_D;

        for(int trg_idx=trg_idx_head;trg_idx<trg_idx_tail;++trg_idx)
        {
            double p[3]={0,0,0};
            double tx=trg[3*vt*stride            + trg_idx];
            double ty=trg[3*vt*stride +   stride + trg_idx];
            double tz=trg[3*vt*stride + 2*stride + trg_idx];

            residual = size_t(src+ 3*stride*vt)%IDEAL_ALIGNMENT;
            if (residual)
                residual = IDEAL_ALIGNMENT - residual;

            // Handle start data if it is not 16-byte aligned
            size_t s;
            // residual = stride;
            for (s=0; s<residual; s++)
            {
                double dX_reg=src[3*stride*vt+           s]-tx;
                double dY_reg=src[3*stride*vt+  stride + s]-ty;
                double dZ_reg=src[3*stride*vt+2*stride + s]-tz;

                double invR = (dX_reg*dX_reg+dY_reg*dY_reg+dZ_reg*dZ_reg);
                if (invR!=0)
                    invR = 1.0/sqrt(invR);

                double cur_pot_x = den[3*stride*vt +           s];
                double cur_pot_y = den[3*stride*vt +  stride + s];
                double cur_pot_z = den[3*stride*vt +2*stride + s];

                double tmp_scalar = (dX_reg*cur_pot_x + dY_reg*cur_pot_y
                    + dZ_reg*cur_pot_z)*invR*invR;
                
                cur_pot_x += tmp_scalar*dX_reg;
                cur_pot_y += tmp_scalar*dY_reg;
                cur_pot_z += tmp_scalar*dZ_reg;

                p[0] += cur_pot_x*invR;
                p[1] += cur_pot_y*invR;
                p[2] += cur_pot_z*invR;
            }

            __m128d txi = _mm_load1_pd (&tx);
            __m128d tyi = _mm_load1_pd (&ty);
            __m128d tzi = _mm_load1_pd (&tz);

            // for (int s=0; s<stride; s++)
            __m128d tempx;
            __m128d tempy;
            __m128d tempz;

            tempx = _mm_setzero_pd();
            tempy = _mm_setzero_pd();
            tempz = _mm_setzero_pd();

            // Load and calculate in groups of SIMD_LEN
            size_t loop_limit = stride-SIMD_LEN_D;
            for (; s <= loop_limit; s += SIMD_LEN_D) {
                __m128d sxj = _mm_load_pd (src+3*stride*vt+         s);
                __m128d syj = _mm_load_pd (src+3*stride*vt+  stride+s);
                __m128d szj = _mm_load_pd (src+3*stride*vt+2*stride+s);

                // this could be vectorized assuming den and q are 16-byte aligned
                __m128d sdenx = _mm_set_pd(
                    den[3*stride*vt +s+1],
                    den[3*stride*vt +s]);

                __m128d sdeny = _mm_set_pd(
                    den[3*stride*vt+stride +s+1],
                    den[3*stride*vt+stride +s]);

                __m128d sdenz = _mm_set_pd (
                    den[3*stride*vt+2*stride +s+1],
                    den[3*stride*vt+2*stride +s]);

                __m128d dX, dY, dZ;
                __m128d dR2;
                __m128d S;

                dX = _mm_sub_pd(txi , sxj);
                dY = _mm_sub_pd(tyi , syj);
                dZ = _mm_sub_pd(tzi , szj);

                sxj = _mm_mul_pd(dX, dX); 
                syj = _mm_mul_pd(dY, dY);
                szj = _mm_mul_pd(dZ, dZ);

                dR2 = _mm_add_pd(sxj, syj);
                dR2 = _mm_add_pd(szj, dR2);

                __m128d zero = _mm_setzero_pd ();
                __m128d is_zero = _mm_cmpeq_pd(dR2, zero);
                
                ///@bug this part of code may be buggy -abtin
                // S = _mm_rsqrt_ps(dR2);
                __m128 dR2_s  = _mm_cvtpd_ps(dR2);
                
                const __m128 approx = _mm_rsqrt_ps( dR2_s );
                const __m128 muls  = _mm_mul_ps(_mm_mul_ps(dR2_s, approx), approx);
                const __m128 three = _mm_set1_ps (3.0f);
                const __m128 half4 = _mm_set1_ps (0.5f);
                __m128 S_s =  _mm_mul_ps(_mm_mul_ps(half4, approx), _mm_sub_ps(three, muls));
                S = _mm_cvtps_pd(S_s);
                S = _mm_andnot_pd (is_zero, S);
                ///@bug

                __m128d dotx = _mm_mul_pd (dX, sdenx);
                __m128d doty = _mm_mul_pd (dY, sdeny);
                __m128d dotz = _mm_mul_pd (dZ, sdenz);

                __m128d dot_sum = _mm_add_pd (dotx, doty);
                dot_sum = _mm_add_pd (dot_sum, dotz);

                dot_sum = _mm_mul_pd (dot_sum, S);
                dot_sum = _mm_mul_pd (dot_sum, S);

                dotx = _mm_mul_pd (dot_sum, dX);
                doty = _mm_mul_pd (dot_sum, dY);
                dotz = _mm_mul_pd (dot_sum, dZ);

                sdenx = _mm_add_pd (sdenx, dotx);
                sdeny = _mm_add_pd (sdeny, doty);
                sdenz = _mm_add_pd (sdenz, dotz);

                sdenx = _mm_mul_pd (sdenx, S);
                sdeny = _mm_mul_pd (sdeny, S);
                sdenz = _mm_mul_pd (sdenz, S);

                tempx = _mm_add_pd (sdenx, tempx);
                tempy = _mm_add_pd (sdeny, tempy);
                tempz = _mm_add_pd (sdenz, tempz);

            }

            _mm_store_pd(tempvalx, tempx); 
            _mm_store_pd(tempvaly, tempy); 
            _mm_store_pd(tempvalz, tempz); 

            for (size_t k = 0; k < SIMD_LEN_D; k++) {
                p[0] += tempvalx[k];
                p[1] += tempvaly[k];
                p[2] += tempvalz[k];
            }

            if (s!=size_t(stride))
                exit(1);
            
            pot[3*vt*stride +            trg_idx] = p[0] * I_PI;
            pot[3*vt*stride +   stride + trg_idx] = p[1] * I_PI;
            pot[3*vt*stride + 2*stride + trg_idx] = p[2] * I_PI;
        }
    }

    return;
}

//////////////////////////////////////////////////////////////////////////
///All to all ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void StokesAlltoAll(const float *src, const float *den, size_t np, float *pot, void*)
{
    float tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    const float *trg(src), *srcPtr(src), *denPtr(den);

#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc, src, den)   
    for (int trg_idx=0; trg_idx<np; ++trg_idx)
    {
        px = 0; py = 0; pz = 0;
        tx=trg[3*trg_idx]; ty=trg[3*trg_idx+1]; tz=trg[3*trg_idx+2];
        
        src = srcPtr; den = denPtr;
        for (int src_idx=0; src_idx<np; ++src_idx)
        {
            dx=*src++ - tx;
            dy=*src++ - ty;
            dz=*src++ - tz;
                
            invR = dx*dx + dy*dy + dz*dz;
            
            if (invR!=0)
                invR = 1.0/sqrt(invR);
            
            cpx = *den++; cpy = *den++; cpz = *den++; 
            
            cc  = dx*cpx + dy*cpy + dz*cpz;
            cc *= invR; cc *= invR;
            
            cpx += cc*dx;
            cpy += cc*dy;
            cpz += cc*dz;
            
            px += cpx*invR;
            py += cpy*invR;
            pz += cpz*invR;
            
        }
        pot[3*trg_idx  ] = px * I_PI;
        pot[3*trg_idx+1] = py * I_PI;
        pot[3*trg_idx+2] = pz * I_PI;
    }
}

void StokesAlltoAll(const double *src, const double *den, size_t np, double *pot, void*)
{
    double tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
    const double *trg(src), *srcPtr(src), *denPtr(den);

#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc, src, den)   
    for (size_t trg_idx=0; trg_idx<np; ++trg_idx)
    {
        px = 0; py = 0; pz = 0;
        tx=trg[3*trg_idx]; ty=trg[3*trg_idx+1]; tz=trg[3*trg_idx+2];
        
        src = srcPtr; den = denPtr;
        for (int src_idx=0; src_idx<np; ++src_idx)
        {
            dx=*src++ - tx;
            dy=*src++ - ty;
            dz=*src++ - tz;
                
            invR = dx*dx + dy*dy + dz*dz;
            
            if (invR!=0)
                invR = 1.0/sqrt(invR);
            
            cpx = *den++; cpy = *den++; cpz = *den++; 
            
            cc  = dx*cpx + dy*cpy + dz*cpz;
            cc *= invR; cc *= invR;
            
            cpx += cc*dx;
            cpy += cc*dy;
            cpz += cc*dz;
            
            px += cpx*invR;
            py += cpy*invR;
            pz += cpz*invR;
            
        }
        pot[3*trg_idx  ] = px * I_PI;
        pot[3*trg_idx+1] = py * I_PI;
        pot[3*trg_idx+2] = pz * I_PI;
    }
}
