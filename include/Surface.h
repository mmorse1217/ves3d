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
#include "GLIntegrator.h"
#include "OperatorsMats.h"
#include "Streamable.h"

#include <queue>
#include <memory>

template <typename ScalarContainer, typename VectorContainer>
class Surface;

template <typename S, typename V>
std::ostream& operator<<(std::ostream& output, const Surface<S, V> &sur);

template <typename ScalarContainer, typename VectorContainer>
class Surface : public Streamable
{
  public:
    typedef typename ScalarContainer::array_type Arr_t;
    typedef typename ScalarContainer::value_type value_type;
    typedef typename ScalarContainer::device_type device_type;
    typedef VectorContainer Vec_t;
    typedef ScalarContainer Sca_t;

    ///@todo add a default constructor
    Surface(const OperatorsMats<Arr_t> &mats, const Vec_t *x_in,
        int upsample = -1, int filter = -1);

    ~Surface();

    void setPosition(const Vec_t& x_in);
    int getNumberOfSurfaces() const { return(getPosition().getNumSubs()); }
    int getShOrder() const { return(getPosition().getShOrder()); }

    Vec_t& getPositionModifiable();
    const Vec_t& getPosition() const;

    const Vec_t& getNormal() const;
    const Sca_t& getAreaElement() const;
    const Sca_t& getMeanCurv() const;
    const Sca_t& getGaussianCurv() const;

    void grad(const Sca_t &f_in, Vec_t &grad_f_out) const;
    void div(const Vec_t &f_in, Sca_t &div_f_out) const;

    void area(Sca_t &area) const;
    void volume(Sca_t &vol) const;
    void getCenters(Vec_t &centers) const;

    void getSmoothedShapePosition(Vec_t &smthd_pos) const;
    void mapToTangentSpace(Vec_t &vec_fld) const;
    void linearizedMeanCurv(const Vec_t &x_new, Sca_t &h_lin) const;

    // From streamable class --------------------------------------------------
    // ------------------------------------------------------------------------
    virtual Error_t pack(std::ostream &os, Format format) const;
    virtual Error_t unpack(std::istream &is, Format format);

    int upsampleFreq() const {return upsample_freq_;}
    int filterFreq() const {return rep_filter_freq_;}
  private:
    Vec_t x_;
    mutable Vec_t normal_;
    mutable Sca_t w_;
    mutable Sca_t h_;
    mutable Sca_t k_;
    mutable Vec_t cu_;
    mutable Vec_t cv_;

    int upsample_freq_;
    int rep_filter_freq_;

    typedef SHTMats<value_type, device_type> Mats_t;
    SHTrans<Sca_t,Mats_t> sht_;
    SHTrans<Sca_t,Mats_t> sht_rep_filter_;
    SHTrans<Sca_t,Mats_t> sht_rep_upsample_;

    GaussLegendreIntegrator<Sca_t> integrator_;

    mutable bool containers_are_stale_;
    mutable bool first_forms_are_stale_;
    mutable bool second_forms_are_stale_;

    void updateFirstForms() const;
    void updateAll() const;
    void checkContainers() const;

    Surface(Surface const& s_in);
    Surface& operator=(const Surface& rhs);

    ///@todo these can be removed, but updateAll should be rewritten
    mutable Sca_t E, F, G;

    mutable std::queue<Sca_t*> scalar_work_q_;
    std::auto_ptr<Sca_t> checkoutSca() const;
    void recycle(std::auto_ptr<Sca_t> scp) const;
    mutable int checked_out_work_sca_;

    mutable std::queue<Vec_t*> vector_work_q_;
    std::auto_ptr<Vec_t> checkoutVec() const;
    void recycle(std::auto_ptr<Vec_t> vcp) const;
    mutable int checked_out_work_vec_;

    void purgeTheWorkSpace() const;

    friend std::ostream& operator<< <Sca_t,Vec_t>(std::ostream& output, const Surface &sur);
};

#include "Surface.cc"

#endif //_SURFACE_H_
