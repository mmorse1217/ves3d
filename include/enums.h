#ifndef _ENUMS_H_
#define _ENUMS_H_

#include <ostream>

///The enum types for the types of device, the insertion operator <<
///is overloaded for this type.
enum DeviceType {CPU, GPU};

///The enum types for the memory copying action, the insertion
///operator << is overloaded for this type.
enum MemcpyKind {MemcpyHostToHost, MemcpyHostToDevice, 
                 MemcpyDeviceToHost, MemcpyDeviceToDevice};

///The enum types for the device related errors, the insertion
///operator << is overloaded for this type.
enum DeviceError {Success, InvalidDevice, SetOnActiveDevice};

///The enum type for the reordering of the points
enum CoordinateOrder {PointMajor, AxisMajor};

enum SolverScheme {Explicit, SemiImplicit};

std::ostream& operator<<(std::ostream& output, const enum DeviceType &DT);
std::ostream& operator<<(std::ostream& output, const enum MemcpyKind &MK);
std::ostream& operator<<(std::ostream& output, const enum DeviceError &err);
std::ostream& operator<<(std::ostream& output, const enum SolverScheme &SS);

#endif //_ENUMS_H_
