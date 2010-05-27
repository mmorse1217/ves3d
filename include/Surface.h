/**
 * @file   Surface.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   Tue Feb  2 13:23:11 2010
 * 
 * @brief  The declaration for the Surface class.
 */

#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "HelperFuns.h"
#include "SHTrans.h"

template <typename ScalarContainer, typename VectorContainer> 
class Surface
{
  public:
    typedef typename ScalarContainer::value_type value_type;
    typedef VectorContainer Vec;
    typedef ScalarContainer Sca;

    ///@todo add a default constructor
    Surface(const Vec& x_in);
    
    void setPosition(const Vec& x_in);
    Vec& getPosition();
    const Vec& getPosition() const;
    const Vec& getNormal() const;
    const Sca& getAreaElement() const;
    const Sca& getMeanCurv() const;
    const Sca& getGaussianCurv() const;

    void updateFirstForms() const;
    void updateAll() const;

    void grad(const Sca &f_in, Vec &grad_f_out) const;
    void div(const Vec &f_in, Sca &div_f_out) const;
 
    void resize(int n_surfs_in);

    void area(Sca &area) const;
    void volume(Sca &vol) const;
    void populate(const Sca &centers);
    void getCenters(Sca &centers) const;

  private:
    VectorContainer x_;
    mutable VectorContainer normal_;

    mutable ScalarContainer w_;
    mutable ScalarContainer h_;
    mutable ScalarContainer k_;

    mutable VectorContainer cu_;
    mutable VectorContainer cv_;
 
    SHTrans<value_type, CPU> sht_;

    mutable bool position_has_changed_outside_;
    mutable bool first_forms_are_stale_;
    mutable bool second_forms_are_stale_;

    void checkContainers() const;
    Surface(Surface<Sca, Vec> const& s_in);
    Surface<Sca, Vec>& operator=(const Surface<Sca, Vec>& rhs);

    ///@todo remove the work vectors
    mutable VectorContainer shc, work_arr;    
    mutable ScalarContainer S1, S2, S3, S4, S5, S6;
    mutable VectorContainer V1, V2;
};

#include "Surface.cc"

#endif //_SURFACE_H_
