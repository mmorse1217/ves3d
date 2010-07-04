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

enum SolverScheme {Explicit, SemiImplicit};

///////////////////////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& output, const enum DeviceType &DT)
{    
    switch (DT)
    {
        case CPU:
            output<<"CPU";
            break;
        case GPU:
            output<<"GPU";
            break;
    }
    
    return output;
}    

std::ostream& operator<<(std::ostream& output, const enum MemcpyKind &MK)
{
    switch (MK)
    {
        case MemcpyHostToHost:
            output<<"MemcpyHostToHost";
            break;

        case MemcpyHostToDevice:
            output<<"MemcpyHostToDevice";
            break;

        case MemcpyDeviceToHost:
            output<<"MemcpyDeviceToHost";
            break;

        case MemcpyDeviceToDevice:
            output<<"MemcpyDeviceToDevice";
            break;
    }
    
    return output;
}    

std::ostream& operator<<(std::ostream& output, const enum DeviceError &err)
{
    switch (err)
    {
        case Success:
            output<<"DeviceError: Success";
            break;
        case InvalidDevice:
            output<<"DeviceError: InvalidDevice";
            break;
        case SetOnActiveDevice:
            output<<"DeviceError: SetOnActiveDevice";
            break;
    }
    
    return output;
}    

std::ostream& operator<<(std::ostream& output, const enum SolverScheme &SS)
{    
    switch (SS)
    {
        case Explicit:
            output<<"Explicit";
            break;
        case SemiImplicit:
            output<<"SemiImplicit";
            break;
    }
    
    return output;
}    

#endif //_ENUMS_H_
