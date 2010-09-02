#include "enums.h"

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

std::ostream& operator<<(std::ostream& output, const enum MonitorReturn &MR)
{    
    switch (MR)
    {
        case StatusOK:
            output<<"StatusOK";
            break;
        case TimeHorizonReached:
            output<<"TimeHorizonReached";
            break;
        case AreaErrorLarge:
            output<<"AreaErrorLarge";
            break;
    }
    
    return output;
}    
