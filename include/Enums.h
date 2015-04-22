/**
 * @file   Enums.h
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   2014-08-26 18:17
 *
 * @brief Enums and macros used in the code and overloaded functions
 * for the enums.
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

#ifndef _ENUMS_H_
#define _ENUMS_H_

#include <iostream>
#include <string>
#include <utility>  // for pair

#define DIM 3

//! typed pi
template <typename T>
inline T PI64(){
    //64   -.------1-------2-------3-------4-------5-------6-------7-------8
    return 3.141592653589793238462643383279502884197169399375105820974944592;
}

///The enum type for the reordering of the points
///PointMajor means [x_1,y_1,z_1,...] ordering
///AxisMajor means [x_1,x_2,...,x_N,...,z_1,...,z_N] ordering
enum CoordinateOrder {PointMajor, AxisMajor, UnknownOrder};

///The linear solver scheme for the vesicle evolution equation
enum SolverScheme {JacobiBlockExplicit,     /* Jacobi iteration + only block tension solve, explicit position */
		   JacobiBlockGaussSeidel,  /* Jacobi iteration + tension solve + position solve              */
		   JacobiBlockImplicit,     /* Jacobi iteration + block coupled solve                         */
		   GloballyImplicit,        /* Fully implicit                                                 */
		   UnknownScheme};          /* Used to signal parsing errors                                  */

///The linear solver scheme for the vesicle evolution equation
enum PrecondScheme {DiagonalSpectral,       /* Only the self preconditioner; diagonal in SH basis */
		    NoPrecond,              /* No preconditioner at all                           */
		    UnknownPrecond};        /* Used to signal parsing errors                      */

///The types of background flow that are supported
enum BgFlowType {ShearFlow,
		 ExtensionalFlow,
		 ParabolicFlow,
		 PeriodicFlow,
		 UserDefinedFlow,        /* Mostly for testing purposes   */
		 UnknownFlow};           /* Used to signal parsing errors */

///DirectEagerEval gives the rotation code the freedom to precompute
///and cache some of the expected results
enum SingularStokesRot {Direct, ViaSpHarm, DirectEagerEval};

///String to enums functions
enum CoordinateOrder EnumifyCoordinateOrder(const char * co);
enum SolverScheme EnumifyScheme(const char * name);
enum PrecondScheme EnumifyPrecond(const char * name);
enum BgFlowType EnumifyBgFlow(const char * name);
enum SingularStokesRot EnumifyStokesRot(const char * name);

std::ostream& operator<<(
    std::ostream& output,
    const enum CoordinateOrder &O);

std::ostream& operator<<(
    std::ostream& output,
    const enum SolverScheme &SS);

std::ostream& operator<<(
    std::ostream& output,
    const enum PrecondScheme &PR);

std::ostream& operator<<(
    std::ostream& output,
    const enum BgFlowType &BG);

std::ostream& operator<<(
    std::ostream& output,
    const enum SingularStokesRot &SS);

#endif //_ENUMS_H_
