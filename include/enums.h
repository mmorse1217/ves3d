#ifndef _ENUMS_H_
#define _ENUMS_H_

#include <ostream>
#include <stack>
#include "Logger.h"

///The enum types for the types of device, the insertion operator <<
///is overloaded for this type.
enum DeviceType {CPU, GPU};

///The enum types for the memory copying action, the insertion
///operator << is overloaded for this type.
enum MemcpyKind {MemcpyHostToHost, MemcpyHostToDevice, 
                 MemcpyDeviceToHost, MemcpyDeviceToDevice};

///The enum type for the reordering of the points
enum CoordinateOrder {PointMajor, AxisMajor};

enum SolverScheme {Explicit, SemiImplicit};

///DirectEagerEval gives the rotation code the freedom to precompute
///and cache some of the expected results
enum SingularStokesRot {Direct, ViaSpHarm, DirectEagerEval}; 

///The enum types for the errors in the Ves3D set of function, the
///insertion operator << is overloaded for this type.

///Errors in the Ves3D code
enum Ves3DErrors {Success,
		  InvalidParameter,
                  UnknownError,
                  InvalidDevice, SetOnActiveDevice, 
                  SolverDiverged, 
                  NoInteraction, InteractionFailed, 
                  RepartitioningFailed,
                  AccuracyError};

typedef enum Ves3DErrors Error_t;

///String to enums
enum SolverScheme EnumifyScheme(const char * name);
enum SingularStokesRot EnumifyStokesRot(const char * name);

///Overloaded insertion operators
std::ostream& operator<<(std::ostream& output, const enum DeviceType &DT);
std::ostream& operator<<(std::ostream& output, const enum MemcpyKind &MK);
std::ostream& operator<<(std::ostream& output, const enum SolverScheme &MR);
std::ostream& operator<<(std::ostream& output, const enum SingularStokesRot &SS);
std::ostream& operator<<(std::ostream& output, const enum Ves3DErrors &err);

//Error handling
class ErrorEvent
{
  public:
    ErrorEvent(const Error_t &err, const char* fun, const char* file, int line);
    
    Error_t err_;
    string funname_;
    string filename_;
    int linenumber_;
};

std::ostream& operator<<(std::ostream& output, const ErrorEvent &ee);

class ErrorHandler
{
  public:
    typedef Error_t (*ErrorCallBack)(const Error_t &);
    
    static ErrorCallBack setErrorCallBack(ErrorCallBack call_back);
    static ErrorEvent submitError(ErrorEvent ee, ErrorCallBack call_back_in = NULL);
    static ErrorEvent submitError(const Error_t &err, const char* fun, 
        const char* file, int line, ErrorCallBack call_back_in = NULL);

    static bool errorStatus();
    static ErrorEvent peekLastError();
    static void popLastError();
    static void clearErrorHist();
    static void printErrorLog();

  private:
    static Error_t ringTheCallBack(ErrorEvent &ee, ErrorCallBack cb = NULL);

    static ErrorCallBack call_back_;
    static stack<ErrorEvent> ErrorStack_;
};

#define QC(err) ( ErrorHandler::submitError(err, __FUNCTION__, __FILE__, __LINE__) )
#define ERRORSTATUS() ( ErrorHandler::errorStatus() )
#define BREAKONERROR() ( {if( !ErrorHandler::errorStatus() ), break;} )
#define PRINTERRORLOG() (ErrorHandler::printErrorLog() )
#define CLEARERRORHIST() (ErrorHandler::clearErrorHist() )

#endif //_ENUMS_H_
