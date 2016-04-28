#include "Enums.h"

enum CoordinateOrder EnumifyCoordinateOrder(const char * co)
{
  std::string ns(co);

  if ( ns.compare(0,4,"Poin") == 0 )
      return PointMajor;
  else if ( ns.compare(0,4,"Axis") == 0 )
      return AxisMajor;
  else
      return UnknownOrder;
}

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
      return UnknownScheme;
}

enum PrecondScheme EnumifyPrecond(const char * name)
{
  std::string ns(name);

  if      ( ns.compare(0,8,"Diagonal") == 0 )
      return DiagonalSpectral;
  else if ( ns.compare(0,9,"NoPrecond") == 0 )
      return NoPrecond;
  else
      return UnknownPrecond;
}

enum BgFlowType EnumifyBgFlow(const char * name)
{
  std::string ns(name);

  if      ( ns.compare(0,5,"Shear") == 0 )
      return ShearFlow;
  else if ( ns.compare(0,9,"Extension") == 0 )
      return ExtensionalFlow;
  else if ( ns.compare(0,9,"Parabolic") == 0 )
      return ParabolicFlow;
  else if ( ns.compare(0,6,"Taylor") == 0 )
      return TaylorVortexFlow;
  else if ( ns.compare(0,6,"Period") == 0 )
      return PeriodicFlow;
  else if ( ns.compare(0,7,"Twister") == 0 )
    return TwisterFlow;
  else if ( ns.compare(0,4,"User") == 0 )
      return UserDefinedFlow;
  else
      return UnknownFlow;
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

enum ReparamType EnumifyReparam(const char * name)
{
  std::string ns(name);
  if ( ns.compare(0,3,"Box") == 0 )
    return BoxReparam;
  else if ( ns.compare(0,5,"PolyK") == 0 )
    return PolyKReparam;
  else
    return UnknownReparam;
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
        case UnknownOrder:
            output<<"UnknownOrder";
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
            output<<"UnknownScheme";
            break;
    }

    return output;
}

std::ostream& operator<<(std::ostream& output, const enum PrecondScheme &PS)
{
    switch (PS)
    {
	case DiagonalSpectral:
            output<<"DiagonalSpectral";
            break;
	case NoPrecond:
            output<<"NoPrecond";
            break;
        default:
            output<<"UnknownPrecond";
            break;
    }

    return output;
}

std::ostream& operator<<(std::ostream& output, const enum BgFlowType &BG)
{
    switch (BG)
    {
	case ShearFlow:
	    output<<"ShearFlow";
	    break;
	case ExtensionalFlow:
	    output<<"ExtensionalFlow";
	    break;
	case ParabolicFlow:
	    output<<"ParabolicFlow";
	    break;
	case TaylorVortexFlow:
	    output<<"TaylorFlow";
	    break;
	case PeriodicFlow:
	    output<<"PeriodicFlow";
	    break;
	case TwisterFlow:
	    output<<"TwisterFlow";
	    break;
	case UserDefinedFlow:
	    output<<"UserDefinedFlow";
	    break;
	default:
	    output<<"UnknownFlow";
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

std::ostream& operator<<(std::ostream& output, const enum ReparamType &RT)
{
    switch (RT)
    {
        case BoxReparam:
            output<<"BoxReparam";
            break;
        case PolyKReparam:
            output<<"PolyKReparam";
            break;
        default:
            output<<"UnknownReparam";
    }

    return output;
}
