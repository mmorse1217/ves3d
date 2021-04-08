#ifndef _SPHERICAL_HARMONICS_H_
#define _SPHERICAL_HARMONICS_H_

#include <matrix.hpp>
#define SHMAXDEG 256

template <class Real>
class SphericalHarmonics{

  public:

    static void SHC2Grid(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>& X, pvfmm::Vector<Real>* X_theta=NULL, pvfmm::Vector<Real>* X_phi=NULL);
    static void SHC2GridPerVesicle(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>& X, pvfmm::Vector<Real>* X_theta=NULL, pvfmm::Vector<Real>* X_phi=NULL);

    static void Grid2SHC(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S);
    static void Grid2SHCPerVesicle(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S);

    static void SHC2GridTranspose(const pvfmm::Vector<Real>& X, long p0, long p1, pvfmm::Vector<Real>& S);

    static void SHC2Pole(const pvfmm::Vector<Real>& S, long p0, pvfmm::Vector<Real>& P);

    static void RotateAll(const pvfmm::Vector<Real>& S, long p0, long dof, pvfmm::Vector<Real>& S_);

    static void RotateTranspose(const pvfmm::Vector<Real>& S_, long p0, long dof, pvfmm::Vector<Real>& S);

    static pvfmm::Vector<Real>& LegendreNodes(long p1);

    static pvfmm::Vector<Real>& LegendreWeights(long p1);

    static pvfmm::Vector<Real>& SingularWeights(long p1);

    static pvfmm::Matrix<Real>& MatFourier(long p0, long p1);

    static pvfmm::Matrix<Real>& MatFourierInv(long p0, long p1);

    static pvfmm::Matrix<Real>& MatFourierGrad(long p0, long p1);

    static std::vector<pvfmm::Matrix<Real> >& MatLegendre(long p0, long p1);

    static std::vector<pvfmm::Matrix<Real> >& MatLegendreInv(long p0, long p1);

    static std::vector<pvfmm::Matrix<Real> >& MatLegendreGrad(long p0, long p1);

    static std::vector<pvfmm::Matrix<Real> >& MatRotate(long p0);

    static void StokesSingularInteg(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>* SLMatrix=NULL, pvfmm::Vector<Real>* DLMatrix=NULL);
    
    static void LaplaceSingularInteg(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>* SLMatrix=NULL, pvfmm::Vector<Real>* DLMatrix=NULL);
    
    static void LaplaceSingularIntegTrans(const pvfmm::Vector<Real>& S, long p0, long p1, pvfmm::Vector<Real>* SLMatrix=NULL, pvfmm::Vector<Real>* DLMatrix=NULL);

    static void SetB0B1();

  private:

    /**
     * \brief Computes all the Associated Legendre Polynomials (normalized) upto the specified degree.
     * \param[in] degree The degree upto which the legendre polynomials have to be computed.
     * \param[in] X The input values for which the polynomials have to be computed.
     * \param[in] N The number of input points.
     * \param[out] poly_val The output array of size (degree+1)*(degree+2)*N/2 containing the computed polynomial values.
     * The output values are in the order:
     * P(n,m)[i] => {P(0,0)[0], P(0,0)[1], ..., P(0,0)[N-1], P(1,0)[0], ..., P(1,0)[N-1],
     * P(2,0)[0], ..., P(degree,0)[N-1], P(1,1)[0], ...,P(2,1)[0], ..., P(degree,degree)[N-1]}
     */
    static void LegPoly(Real* poly_val, const Real* X, long N, long degree);

    static void LegPolyDeriv(Real* poly_val, const Real* X, long N, long degree);

    template <bool SLayer, bool DLayer>
    static void StokesSingularInteg_(const pvfmm::Vector<Real>& X0, long p0, long p1, pvfmm::Vector<Real>& SL, pvfmm::Vector<Real>& DL);

    template <bool SLayer, bool DLayer>
    static void LaplaceSingularInteg_(const pvfmm::Vector<Real>& X0, long p0, long p1, pvfmm::Vector<Real>& SL, pvfmm::Vector<Real>& DL);

    template <bool SLayer=true, bool DLayer=false>
    static void LaplaceSingularIntegTrans_(const pvfmm::Vector<Real>& X0, long p0, long p1, pvfmm::Vector<Real>& SL, pvfmm::Vector<Real>& DL);

    static struct MatrixStorage{
      MatrixStorage(int size){
        Qx_ .resize(size);
        Qw_ .resize(size);
        Sw_ .resize(size);
        Mf_ .resize(size*size);
        Mdf_.resize(size*size);
        Ml_ .resize(size*size);
        Mdl_.resize(size*size);
        Mr_ .resize(size);
        Mfinv_ .resize(size*size);
        Mlinv_ .resize(size*size);
        B0s.resize(omp_get_max_threads());
        B1s.resize(omp_get_max_threads());
      }
      std::vector<pvfmm::Vector<Real> > Qx_;
      std::vector<pvfmm::Vector<Real> > Qw_;
      std::vector<pvfmm::Vector<Real> > Sw_;
      std::vector<pvfmm::Matrix<Real> > Mf_ ;
      std::vector<pvfmm::Matrix<Real> > Mdf_;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Ml_ ;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Mdl_;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Mr_;
      std::vector<pvfmm::Matrix<Real> > Mfinv_ ;
      std::vector<std::vector<pvfmm::Matrix<Real> > > Mlinv_ ;
      std::vector< pvfmm::Vector<Real> > B0s, B1s;
    } matrix;

};

template<> SphericalHarmonics<double>::MatrixStorage SphericalHarmonics<double>::matrix(SHMAXDEG);
template<> SphericalHarmonics<float >::MatrixStorage SphericalHarmonics<float>::matrix(SHMAXDEG);

#include "SphericalHarmonics.cc"

#endif // _SPHERICAL_HARMONICS_H_

