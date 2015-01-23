/**
 * @file
 * @author Rahimian, Abtin <arahimian@acm.org>
 * @date   2014-08-26 17:49
 *
 * @brief Error handling class (singleton) and macros
 */

/*
 * Copyright (c) 2014, Abtin Rahimian
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _ERROR_H_
#define _ERROR_H_

#include <iostream>
#include <stack>
#include <string>

class ErrorEvent
{
  public:
    //! Errors in the Ves3D code
    enum Error_t {Success		 = 0,
		  UnknownError		 = 1,
		  ArithmeticError	 = 2,
		  AssertionError	 = 3,
		  EnvironmentError	 = 4,	/* IO, file, etc */
		  MemoryError		 = 5,
		  ReferenceError	 = 6,
		  NotImplementedError	 = 7,
		  //
                  InvalidParameterError  = 10,
		  //
                  InvalidDeviceError	 = 20,
		  SetOnActiveDeviceError = 21,
		  //
                  AccuracyError		 = 30,
                  DivergenceError	 = 31,
		  //
                  InteractionError	 = 40,
		  RepartitioningError	 = 41};

    //! lightweight error event that can be handled by
    //! ErrorHandler
    ErrorEvent(const Error_t &err, const char* fun,
        const char* file, int line);

    Error_t      err_;
    std::string  funname_;
    std::string  filename_;
    int          linenumber_;

    operator bool(){return bool(err_);}
};

typedef ErrorEvent::Error_t Error_t;

//! Singleton class for error handling
class ErrorHandler
{
  public:
    //! Callback function pointer type
    typedef Error_t (*ErrorCallBack)(const Error_t &);

    //! Updates the currecnt callback and returns the old one
    static ErrorCallBack setErrorCallBack(ErrorCallBack call_back);

    //! Raise an error with optional callback
    static ErrorEvent submitError( ErrorEvent ee,
        ErrorCallBack call_back_in = NULL);
    static ErrorEvent submitError(const Error_t &err, const char* fun,
        const char* file, int line, ErrorCallBack call_back_in = NULL);

    // Error stack lookup
    static bool errorStatus();
    static ErrorEvent peekLastError();
    static void popLastError();
    static void clearErrorHist();
    static void printErrorLog();

  private:
    static Error_t ringTheCallBack(ErrorEvent &ee,
        ErrorCallBack cb = NULL);

    static ErrorCallBack call_back_;
    static std::stack<ErrorEvent> ErrorStack_;
};

std::ostream& operator<<(std::ostream& output,
    const ErrorEvent &ee);

std::ostream& operator<<(std::ostream& output,
    const Error_t &err_t);

// Error macros
//! CHK(err) with err as Error_t enum, raises an error if err is not Success
// defineing chk_expr to avoid multiple evaluation of expr
static Error_t chk_expr(ErrorEvent::UnknownError);

#define CHK(expr) (                            \
		   ( chk_expr = expr ),	       \
		   chk_expr &&		       \
        ErrorHandler::submitError(chk_expr,    \
            __FUNCTION__,                      \
            __FILE__,                          \
            __LINE__                           \
                                  ))

#define CHK_CB(expr,callback) (						\
			       ( chk_expr = expr ),			\
			       chk_expr &&				\
			       ErrorHandler::submitError(chk_expr,	\
							 __FUNCTION__,	\
							 __FILE__,	\
							 __LINE__,	\
							 callback	\
							 ))

#define SET_ERR_CALLBACK(cb) ( ErrorHandler::setErrorCallBack(cb) )
#define ERRORSTATUS() ( ErrorHandler::errorStatus() )
#define BREAKONERROR() ( {if( !ErrorHandler::errorStatus() ) break;} )
#define PRINTERRORLOG() (ErrorHandler::printErrorLog() )
#define CLEARERRORHIST() (ErrorHandler::clearErrorHist() )

#endif //_ERROR_H_
