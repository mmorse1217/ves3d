#ifndef _EVALVELOCITY_H_
#define _EVALVELOCITY_H_

#include "Enums.h"
#include "Surface.h"
#include "InterfacialForce.h"
#include "BgFlowBase.h"

template<typename Scalar,
         typename Vector,
         typename StokesEvaluator>
class EvalVelocity
{
  private:
    typedef Surface<Scalar, Vector> Sur_t;
    typedef typename Scalar::value_type value_type;


  public:
    EvalVelocity(const StokesEvaluator &stokes,
        const BgFlowBase<Vector> &vInf,
        OperatorsMats<Scalar> &mats,
        value_type bending_modulus);
    ~EvalVelocity();

    Error_t operator()(const Vector &x_src,
        const Scalar &tension, const Vector &x_eval,
        Vector &vel);

  private:
    const StokesEvaluator &stokes_;
    const BgFlowBase<Vector> &vInf_;
    OperatorsMats<Scalar> &mats_;
    Sur_t *S_ptr_; //pointer, to avoid construction of the surface in
                   //the constructor
    InterfacialForce<Sur_t> Force_;
    Scalar quad_weights_;

    //Work space
    Vector Fb, Fs;
    Vector all_src, all_den, all_pot;
};

#include "EvalVelocity.cc"

#endif //_EVALVELOCITY_H_
