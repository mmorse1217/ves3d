
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
	case ErrorEvent::UnknownError:
            output<<"UnknownError";
            break;
	case ErrorEvent::ArithmeticError:
            output<<"ArithmeticError";
            break;
	case ErrorEvent::AssertionError:
            output<<"AssertionError";
            break;
	case ErrorEvent::EnvironmentError:
            output<<"EnvironmentError";
            break;
	case ErrorEvent::MemoryError:
            output<<"MemoryError";
            break;
	case ErrorEvent::ReferenceError:
            output<<"ReferenceError";
            break;
	case ErrorEvent::NotImplementedError:
            output<<"NotImplementedError";
            break;
	case ErrorEvent::InvalidParameterError:
            output<<"InvalidParameterError";
            break;
	case ErrorEvent::InvalidDeviceError:
            output<<"InvalidDeviceError";
            break;
	case ErrorEvent::SetOnActiveDeviceError:
            output<<"SetOnActiveDeviceError";
            break;
	case ErrorEvent::AccuracyError:
            output<<"AccuracyError";
            break;
	case ErrorEvent::DivergenceError:
            output<<"DivergenceError";
            break;
	case ErrorEvent::InteractionError:
            output<<"InteractionError";
            break;
	case ErrorEvent::RepartitioningError:
            output<<"RepartitioningError";
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
        return ( cb(ee) );
    else if ( ErrorHandler::call_back_ != NULL )
        return ( ErrorHandler::call_back_(ee) );
    else
    {
        //not using CERR_LOC b/c of duplicate file info
        CERR("No callback is set to handle error event:\n"<<ee);
    }

    return ErrorEvent::UnknownError;
}
