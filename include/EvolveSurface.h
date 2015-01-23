#ifndef _EVOLVESURFACE_H_
#define _EVOLVESURFACE_H_

#include "Logger.h"
#include "Enums.h"
#include "Error.h"

#include "Device.h"
#include "Vectors.h"
#include "Surface.h"
#include "Surface.h"
#include "OperatorsMats.h"
#include "Parameters.h"
#include "BgFlowBase.h"
#include "InterfacialVelocity.h"

//The default arguments classes for the template
#include "BgFlow.h"
#include "VesInteraction.h"
#include "Repartition.h"
#include "Monitor.h"

/**
 * EvolveSurface uses a simple Euler time stepping (explicit,
 * implicit) method to update the surface. Major components of
 * simulation (BgFlow, Interaction, Repartition, etc.) are arguments
 * to this.
 *
 * The pseudocode of the Evolve method is:
 *  \code
 *  while ( t < time_horizon && [no_error])
 *  {
 *     . update surfaces
 *     . reparametrize
 *     . repartition
 *     . monitor
 *  }
 * \endcode
 */
template<typename T,
         typename DT,
         const DT &DEVICE,
         typename Interact = VesInteraction<T>,
         typename Repart   = Repartition<T> >
class EvolveSurface
{
  public:
    typedef T value_type;
    typedef DT device_type;
    typedef Scalars<T, DT, DEVICE> Sca_t;
    typedef typename Sca_t::array_type Arr_t;
    typedef Vectors<T, DT, DEVICE> Vec_t;
    typedef Surface<Sca_t, Vec_t> Sur_t;

    typedef BgFlowBase<Vec_t> BgFlow_t;
    typedef MonitorBase<EvolveSurface> Monitor_t;
    typedef Interact Interaction_t;
    typedef Repart Repartition_t;

    typedef OperatorsMats<Arr_t> Mats_t;
    typedef Parameters<T> Params_t;

    ///external function pointer types
    typedef typename Interact::InteractionFun_t InteractionFun_t;
    typedef typename Repart::GlobalRepart_t GlobalRepart_t;

    typedef InterfacialVelocity<Sur_t, Interaction_t> IntVel_t;
    typedef Error_t (IntVel_t::*Scheme_t)(const value_type &);

    EvolveSurface(Params_t &params, Mats_t &mats, Vec_t &x0, BgFlow_t *vInf,
        Monitor_t *M = NULL, Interaction_t *I = NULL, Repartition_t *R=NULL
        );

    ~EvolveSurface();

    Error_t Evolve();

    Params_t &params_;
    Mats_t &mats_;

    Sur_t *S_;
    BgFlow_t *vInf_;
    Monitor_t *monitor_;
    Interaction_t *interaction_;
    Repartition_t *repartition_;
    IntVel_t *F_;

    bool objsOnHeap_[3];//for monitor, interaction, and repartition
};

#include "EvolveSurface.cc"

#endif //_EVOLVESURFACE_H_
