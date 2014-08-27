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
#include <stack>
#include <string>

using namespace std;

#define DIM 3

///The enum types for the types of device, the insertion operator <<
///is overloaded for this type.
enum DeviceType {CPU, GPU};

///The enum types for the memory copying action, the insertion
///operator << is overloaded for this type.
enum MemcpyKind {MemcpyHostToHost, MemcpyHostToDevice,
                 MemcpyDeviceToHost, MemcpyDeviceToDevice};

///The enum type for the reordering of the points
///PointMajor means [x_1,y_1,z_1,...] ordering
///AxisMajor means [x_1,x_2,...,x_N,...,z_1,...,z_N] ordering
enum CoordinateOrder {PointMajor, AxisMajor};

///The linear solver scheme for the vesicle evolution equation
enum SolverScheme {Explicit, BlockImplicit, GloballyImplicit};

///DirectEagerEval gives the rotation code the freedom to precompute
///and cache some of the expected results
enum SingularStokesRot {Direct, ViaSpHarm, DirectEagerEval};

///The enum types for the errors in the Ves3D set of function, the
///insertion operator << is overloaded for this type.

///Errors in the Ves3D code
enum Ves3DErrors {Success=0,
                  InvalidParameter,
                  UnknownError,
                  InvalidDevice, SetOnActiveDevice,    //Device error
                  SolverDiverged,
                  NoInteraction, InteractionFailed,
                  NoRepartition, RepartitioningFailed,
                  AccuracyError};

typedef enum Ves3DErrors Error_t;

///String to enums functions
enum SolverScheme EnumifyScheme(const char * name);
enum SingularStokesRot EnumifyStokesRot(const char * name);

///Overloaded insertion operators
std::ostream& operator<<(
    std::ostream& output,
    const enum DeviceType &DT);

std::ostream& operator<<(
    std::ostream& output,
    const enum MemcpyKind &MK);

std::ostream& operator<<(
    std::ostream& output,
    const enum CoordinateOrder &O);

std::ostream& operator<<(
    std::ostream& output,
    const enum SolverScheme &MR);

std::ostream& operator<<(
    std::ostream& output,
    const enum SingularStokesRot &SS);

std::ostream& operator<<(
    std::ostream& output,
    const enum Ves3DErrors &err);

#endif //_ENUMS_H_
