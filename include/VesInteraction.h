#ifndef _VESINTERACTION_H_
#define _VESINTERACTION_H_

#include "HelperFuns.h"

#define I_PI 1.0/M_PI/8.0

template<typename VecContainer>
class VesInteraction
{
  public:
    void operator()(const VecContainer &position, VecContainer &density, 
        VecContainer &potential) const 
    {
        size_t np(position.getNumSubs() * position.getStride());
        typename VecContainer::value_type tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc;
        
        typename VecContainer::value_type const *src(position.begin());
        typename VecContainer::value_type const *trg(position.begin());
        typename VecContainer::value_type *den(density.begin());
        typename VecContainer::value_type *pot(potential.begin());

#pragma omp parallel for private(tx, ty, tz, px, py, pz, dx, dy, dz, invR, cpx, cpy, cpz, cc)
        for (size_t trg_idx=0;trg_idx<np; ++trg_idx)
        {
            px = 0; py = 0; pz = 0;
            
            tx= *trg++;
            ty= *trg++;
            tz= *trg++;
            
//             for (size_t src_idx=0; src_idx<np; ++src_idx)
//             {
//                 dx= *(src + 3 * src_idx)-tx;
//                 dy= *(src + 3 * src_idx + 1)-ty;
//                 dz= *(src + 3 * src_idx + 2)-tz;

//                 invR = dx*dx; invR+= dy*dy; invR+= dz*dz;
//                 if (invR!=0)
//                     invR = 1.0/sqrt(invR);
            
//                 cpx = *(den + 3 * src_idx); 
//                 cpy = *(den + 3 * src_idx + 1);  
//                 cpz = *(den + 3 * src_idx + 2);  
                
//                 cc  = dx*cpx; cc += dy*cpy; cc += dz*cpz;
//                 cc *= invR; cc *= invR;

//                 cpx += cc*dx; cpy += cc*dy; cpz += cc*dz;
//                 px += cpx*invR; py += cpy*invR; pz += cpz*invR;
//             }
            
            *pot++ = px * I_PI;
            *pot++ = py * I_PI;
            *pot++ = pz * I_PI;
        }
    }
};

#endif // _VESINTERACTION_H_
