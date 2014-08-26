#include "Enums.h"

enum SolverScheme EnumifyScheme(const char * name)
{
  string ns(name);

  if ( ns.compare(0,8,"Explicit") == 0 )
    return Explicit;
  else if ( ns.compare(0,5,"Block") == 0 )
    return BlockImplicit;
  else
    return GloballyImplicit;
}
enum SingularStokesRot EnumifyStokesRot(const char * name)
{
  string ns(name);
  if ( ns.compare(0,9,"ViaSpHarm") == 0 )
    return ViaSpHarm;
  else if ( ns.compare(0,15,"DirectEagerEval") == 0 )
    return DirectEagerEval;
  else
    return Direct;
}

//
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

std::ostream& operator<<(std::ostream& output,
    const enum CoordinateOrder &O)
{
    switch (O)
    {
        case PointMajor:
            output<<"PointMajor";
            break;
        case AxisMajor:
            output<<"AxisMajor";
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
        case BlockImplicit:
            output<<"BlockImplicit";
            break;
        case GloballyImplicit:
            output<<"GloballyImplicit";
            break;
    }

    return output;
}

std::ostream& operator<<(std::ostream& output, const enum SingularStokesRot &SS)
{
    switch (SS)
    {
        case Direct:
            output<<"Direct";
            break;
        case ViaSpHarm:
            output<<"ViaSpHarm";
            break;
        case DirectEagerEval:
            output<<"DirectEagerEval";
            break;
    }

    return output;
}

std::ostream& operator<<(std::ostream& output, const enum Ves3DErrors &err)
{
    switch (err)
    {
        case Success:
            output<<"Success";
            break;
        case InvalidParameter:
            output<<"InvalidParameter";
	    break;
        case UnknownError:
            output<<"UnkownError";
            break;
        case InvalidDevice:
            output<<"InvalidDevice";
            break;
        case SetOnActiveDevice:
            output<<"SetOnActiveDevice";
            break;
        case SolverDiverged:
            output<<"SolverDiverged";
            break;
        case NoInteraction:
            output<<"NoInteraction";
            break;
        case InteractionFailed:
            output<<"InteractionFailed";
            break;
        case RepartitioningFailed:
            output<<"RepartitioningFailed";
            break;
        case AccuracyError:
            output<<"AccuracyError";
            break;
        default:
            output<<err
                  <<" [The string for the given enum type is not known"
                  <<", update the insertion operator]";
    }

    return output;
}
