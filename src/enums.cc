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

std::ostream& operator<<(std::ostream& output, const enum Ves3DErrors &err)
{    
    switch (err)
    {
        case Success:
            output<<"Success";
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
            output<<err<<" [The string for the given enum type is not known, update the insertion operator]";
    }
    
    return output;
}


//Error handling
ErrorEvent::ErrorEvent(const Error_t &err, const char* fun, const char* file, int line) :
    err_(err),
    funname_(fun),
    filename_(file),
    linenumber_(line)
{}

std::ostream& operator<<(std::ostream& output, const ErrorEvent &ee)
{
    output<<"\n       Error         : "<<ee.err_
          <<"\n       Function name : "<<ee.funname_
          <<"\n       File name     : "<<ee.filename_
          <<"\n       Line number   : "<<ee.linenumber_<<"\n";
    
    return output;
}

stack<ErrorEvent> ErrorHandler::ErrorStack_;
ErrorHandler::ErrorCallBack ErrorHandler::call_back_(NULL);

ErrorHandler::ErrorCallBack ErrorHandler::setErrorCallBack(ErrorHandler::ErrorCallBack call_back)
{
    ErrorCallBack tmp(call_back_);
    ErrorHandler::call_back_ = call_back;
    return tmp;
}

ErrorEvent ErrorHandler::submitError(ErrorEvent ee, ErrorCallBack call_back_in)
{
    if (ee.err_ != Success )
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
    stack<ErrorEvent> tmp_stack;
    
    COUT(" ===================================="<<endl);
    if ( ErrorStack_.empty() ) 
        COUT("  Error log is clean."<<endl);
    else
        COUT("   Error Log\n ------------------------------------"<<endl);

    while ( !ErrorStack_.empty() )
    {
        tmp_stack.push(ErrorStack_.top());
        ErrorStack_.pop();
    }

    while ( !tmp_stack.empty() )
    { 
        ErrorStack_.push(tmp_stack.top());
        COUT(ErrorStack_.top()<<endl);
        COUT(" ------------------------------------"<<endl);
        tmp_stack.pop();
    }

    COUT(" ===================================="<<endl);
}

Error_t ErrorHandler::ringTheCallBack(ErrorEvent &ee, ErrorHandler::ErrorCallBack cb)
{
    if ( cb != NULL )
        return ( cb(ee.err_) );
    else if ( ErrorHandler::call_back_ != NULL )
        return ( ErrorHandler::call_back_(ee.err_) ); 
    else
    {
        cerr<<"\n  An Error was encountered during the function invocation in:\n"<<ee<<endl;
    }
    
    return UnknownError;
}
