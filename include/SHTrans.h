#ifndef _SHTRANS_H_
#define _SHTRANS_H_

/**
 * Spherical Harmonics Transform (SHT) class. The template parameter
 * <code>Container</code> is assumed to have a static method <code>
 * Container::getDevice()</code> that returns its associated
 * device. The returned device is assumed to have <code> Malloc(),
 * Memcpy(), Free(), gemm(), Transpose(),</code> and <code>
 * ax()</code> methods. See the Device class for further
 * information. The <code>Mats</code> template parameter holds the
 * Analysis and Synthesis matrices of discrete Fourier transform,
 * discrete Legendre transform, and first and second derivatives of
 * latitude(Legendre) variable. See SHTMats for more information.
 */
template<typename Container, typename Mats>
class SHTrans
{
  private:
    typedef typename Container::value_type value_type;

  public:
     SHTrans(int sh_order_in, Mats &mats, int filter_freq = -1);
    ~SHTrans();

    int getShOrder() const {return(p);}

    /**
     * The Synthesis method. The input is in real space and the output
     * is the spherical harmonics coefficients (shc).  The input,
     * work, and shc are all of the same size.
     *
     * The shc has the real spherical harmonic's coefficent:
     *  [Y_0,0 Y_1,0... Y_p,0 Y1,-1 Y_2,-1 ... Y_p,-1 Y_1,1 ... Y_p,1 ...]
     */
    void forward(const Container &in, Container &work, Container &shc) const;

    //Analysis Methods (all are wrappers around the private back method):
    void backward(const Container &shc, Container &work,
        Container &out) const;
    void backward_du(const Container &shc, Container &work,Container &out) const;
    void backward_dv (const Container &shc, Container &work, Container &out) const;
    void backward_d2u(const Container &shc, Container &work, Container &out) const;
    void backward_d2v(const Container &shc, Container &work, Container &out) const;
    void backward_duv(const Container &shc, Container &work, Container &out) const;

    void FirstDerivatives(const Container &in, Container &work,
        Container &shc, Container &du, Container &dv) const;

    void lowPassFilter(const Container &in, Container &work,
        Container &shc, Container &out) const;

    void collectSameOrder(const Container &in, Container &out) const;
    void collectSameFreq(const Container &in, Container &out) const;

    void ScaleFreq(const value_type *shc_in, int n_funs,
        const value_type* scaling_coeff, value_type *shc_out) const;

  private:
    typedef typename Container::device_type device_type;
    const typename Container::device_type &device_;

    Mats &mats_;

    static const value_type alpha_ = (value_type) 1.0;
    static const value_type beta_  = (value_type) 0.0;

    int p;
    int dft_size;

    /**
     * The forward and inverse Legendre transform. The case of forward
     * and backward depends on the matrix (set) <code>trans</code>
     * passed to the method.
     *
     * @param trans The set of transform matrices saves one after
     * another in one long array.
     * @param m The first dimension of the largest trans matrix
     * @param n The second dimension of input/output
     * @param k The second dimension of trans matrix and the first
     * dimension of input
     * @param mf The flag indicating whether the trans matrices size m
     * decrease or not (used when called in the forward transform)
     * @param nf The flag indicating whether the input matrix size
     * decreases (not used)
     * @param kf The flag indicating whether the trans matrices size k
     * decrease or not (used in the call from the inverse transforms)
     */
    void DLT(value_type *trans, const value_type *inputs,
        value_type *outputs, int m, int n , int k, int mf,
        int nf, int kf) const;

    ///The inverse transform.
    void back(const value_type *inputs, value_type *work_arr,
        int n_funs, value_type *outputs, value_type *trans,
        value_type *dft) const;

    value_type* filter_coeff_;
};

#include "SHTrans.cc"
#endif //_SHTRANS_H_
