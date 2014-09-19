/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @revision $Revision$
 * @tags $Tags$
 * @date $Date$
 *
 * @brief The interface to spherical harmonics transform and factory
 * functions to generate transformation class.
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _SPHARM_H_
#define _SPHARM_H_

#include <utility> //pair
#include <cstddef> //size_t
#include "Error.h"
#include "tr1.h"   //shared_ptr


//! Sentinel grid-size for spherical harmonic
#define EMPTY_GRID std::make_pair(-1,-1)

//! Grid function for spherical harmonics (longitude,lattitude) grid
inline std::pair<int, int> SpharmGridDim(int sh_order)
{
    return((sh_order >= 0) ?
        std::make_pair(sh_order+1, 2*sh_order+2) :
        EMPTY_GRID);
}

/**
 * Spherical harmonics interface declaration.
 *
 * @param T content type (float or double)
 * @note  Derived classes should implement static factory and recycle
 * functions (although not explicit in the interface).
 */
template<typename T> //T is content type
class Spharm{
  public:
    typedef T value_type;

    virtual ~Spharm() = 0;

    /**
     * Transform function values in f to coefficients in c.
     *
     * @param[in]  f function samples in real space over the grid
     * @param[out] c spherical harmonic coefficient
     * @param[in]  nf number of functions.
     * @param[in]  sh_order spherical harmonic order of function.
     * @param[in]  stride length of each function (deduced from
     *             sh_order if default)
     * @return     Error_t enum to indicate success or failure
     *
     * @note       The expected length of f and c is stride*nf.
     * @see        Evaluate()
     */
    virtual Error_t Coeffs(const T* f, T* c, size_t nf,
        int sh_order = -1, size_t stride = 0) = 0;

    /**
     * Evaluate the coefficients on the grid and returns function values
     *
     * @param[in]  c spherical harmonics coefficients
     * @param[out] f function evaluated on the grid
     * @return       Error_t for success or failure
     * @see          Coeffs()
     */
    virtual Error_t Evaluate(const T* c, T* f, size_t nf,
        int sh_order = -1, size_t stride = 0) = 0;
};

//some compilers complain about undefined deconstructor
template<typename T>
Spharm<T>::~Spharm(){}

/**
 * Factory function for derived classes. This implies that derived
 * classes have static member Factory().
 *
 * @param[in] sh_order spherical harmonic order the transform is needed
 * @return    pointer to an instance of Derived class
 */
template<typename Derived>
inline tr1::shared_ptr<Spharm<typename Derived::value_type> >SpharmFactory(
    int sh_order){
    return( Derived::Factory(sh_order) );
}

#endif //_SPHARM_H_
