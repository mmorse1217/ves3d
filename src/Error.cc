#include "Error.h"
#include "Logger.h"

/*
 * Error event
 */
ErrorEvent::ErrorEvent(const Error_t &err, const char* fun,
    const char* file, int line) :
    err_(err),
    funname_(fun),
    filename_(file),
    linenumber_(line)
{}

std::ostream& operator<<(std::ostream& output, const ErrorEvent &ee)
{
    output<<"    Error         : "<<ee.err_
          <<"\n    Function name : "<<ee.funname_
          <<"\n    File name     : "<<ee.filename_
          <<"\n    Line number   : "<<ee.linenumber_;

    return output;
}

std::ostream& operator<<(std::ostream& output, const Error_t &err_t)
{
    switch (err_t)
    {
        case ErrorEvent::Success:
            output<<"Success";
            break;
        case ErrorEvent::InvalidParameter:
            output<<"InvalidParameter";
	    break;
        case ErrorEvent::UnknownError:
            output<<"UnkownError";
            break;
        case ErrorEvent::InvalidDevice:
            output<<"InvalidDevice";
            break;
        case ErrorEvent::SetOnActiveDevice:
            output<<"SetOnActiveDevice";
            break;
        case ErrorEvent::SolverDiverged:
            output<<"SolverDiverged";
            break;
        case ErrorEvent::NoInteraction:
            output<<"NoInteraction";
            break;
        case ErrorEvent::InteractionFailed:
            output<<"InteractionFailed";
            break;
        case ErrorEvent::RepartitioningFailed:
            output<<"RepartitioningFailed";
            break;
        case ErrorEvent::AccuracyError:
            output<<"AccuracyError";
            break;
        default:
            output<<err_t
                  <<" [The string for the given enum type is not known"
                  <<", update the insertion operator]";
    }

    return output;
}

/*
 * Error handler
 */
std::stack<ErrorEvent> ErrorHandler::ErrorStack_;
ErrorHandler::ErrorCallBack ErrorHandler::call_back_(NULL);

ErrorHandler::ErrorCallBack ErrorHandler::setErrorCallBack(
    ErrorHandler::ErrorCallBack call_back)
{
    ErrorCallBack tmp(call_back_);
    ErrorHandler::call_back_ = call_back;
    return tmp;
}

ErrorEvent ErrorHandler::submitError(ErrorEvent ee,
    ErrorCallBack call_back_in)
{
    if (ee.err_ != ErrorEvent::Success )
    {
        ErrorHandler::ErrorStack_.push(ee);
        ErrorHandler::ringTheCallBack(ee, call_back_in);
    }
    return ee;
}

ErrorEvent ErrorHandler::submitError(const Error_t &err, const char* fun,
    const char* file, int line, ErrorCallBack call_back_in)
{
    return( ErrorHandler::submitError(ErrorEvent(err, fun, file, line),
                                          call_back_in) );
}

bool ErrorHandler::errorStatus()
{
    return( ErrorStack_.empty() );
}

ErrorEvent ErrorHandler::peekLastError()
{
    return( ErrorStack_.top() );
}

void ErrorHandler::popLastError()
{
    return( ErrorStack_.pop() );
}

void ErrorHandler::clearErrorHist()
{
    while ( !ErrorStack_.empty() )
        ErrorStack_.pop();
}

void ErrorHandler::printErrorLog()
{
    std::stack<ErrorEvent> tmp_stack;

    if ( ErrorStack_.empty() ){
        COUT("=============================================="
            <<"=======================================\n"
            <<"Error log is clean.\n"
            <<"=============================================="
            <<"=======================================");
        return;
    }

    CERR("=============================================="
        <<"=======================================");
    CERR("Error Log");

    while ( !ErrorStack_.empty() )
    {
        tmp_stack.push(ErrorStack_.top());
        ErrorStack_.pop();
    }

    while ( !tmp_stack.empty() )
    {
        ErrorStack_.push(tmp_stack.top());
        CERR("---------------------------------------------"
            <<"----------------------------------------");
        CERR(ErrorStack_.top());
        tmp_stack.pop();
    }

    CERR("============================================"
        <<"=========================================");
}

Error_t ErrorHandler::ringTheCallBack(ErrorEvent &ee,
    ErrorHandler::ErrorCallBack cb)
{
    if ( cb != NULL )
        return ( cb(ee.err_) );
    else if ( ErrorHandler::call_back_ != NULL )
        return ( ErrorHandler::call_back_(ee.err_) );
    else
    {
        //not using CERR b/c of duplicate file info
        CERR_LOC("No callback is set to handle error event:\n"<<ee,"",NULL);
    }

    return ErrorEvent::UnknownError;
}
