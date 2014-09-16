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
