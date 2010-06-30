#ifndef _ENUMS_H_
#define _ENUMS_H_

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

#endif //_ENUMS_H_
