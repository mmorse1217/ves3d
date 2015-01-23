#include "Enums.h"

enum SolverScheme EnumifyScheme(const char * name)
{
  std::string ns(name);

  if      ( ns.compare(0,12,"JacobiBlockE") == 0 )
      return JacobiBlockExplicit;
  else if ( ns.compare(0,12,"JacobiBlockG") == 0 )
      return JacobiBlockGaussSeidel;
  else if ( ns.compare(0,12,"JacobiBlockI") == 0 )
      return JacobiBlockImplicit;
  else if ( ns.compare(0,6 ,"Global") == 0 )
      return GloballyImplicit;
  else
      return UnkownScheme;
}
enum SingularStokesRot EnumifyStokesRot(const char * name)
{
  std::string ns(name);
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
	case JacobiBlockExplicit:
            output<<"JacobiBlockExplicit";
            break;
	case JacobiBlockGaussSeidel:
            output<<"JacobiBlockGaussSeidel";
            break;
        case JacobiBlockImplicit:
            output<<"JacobiBlockImplicit";
            break;
        case GloballyImplicit:
            output<<"GloballyImplicit";
            break;
        default:
            output<<"UnkownScheme";
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
