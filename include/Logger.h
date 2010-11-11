#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream>
//#include <sys/time.h>
#include <map>
#include <stack>
#include <algorithm>
#include <string>
#include <omp.h>

using namespace std;

///A log event recorded by the Logger class.
struct LogEvent
{
    string fun_name;
    unsigned long int num_calls;
    double time;
    double flop;
    double flop_rate;
};

///The print function for the LogEvent class.
template<typename T>
void PrintLogEvent(const pair<T, LogEvent>& ev);

///The enum type for the format of the report generated by
///Logger::Report().
enum ReportFormat {SortFunName, SortNumCalls, SortTime, 
                   SortFlop, SortFlopRate};

class Logger
{
  public:
    ///Returns the current wall-time.
    static double Now();
    
    ///Gets the current wall-time, saves it in a stack and also
    ///returns it as output.
    static double Tic();
    
    ///Gets the current wall-time, pops the corresponding Tic() value
    ///form the stack and returns the difference.
    static double Toc();

    ///Records accumulative log events corresponding to functions
    ///calling this method.
    static void Record(string fun_name, string prefix, double time,
        double flop);
    
    ///Reports all accumulative data corresponding to all the calls to
    ///the Record() function.
    static void Report(enum ReportFormat rf);

    ///The setter function of the log file.
    static void SetLogFile(string file_name);
    
    ///The method to log events, it writes the event to the file
    ///log_file.
    static void Log(const char *event);

    ///Clears the slate of profiler
    static void PurgeProfileHistory();

    static double GetFlops();
  private:
    ///The constructor, since the class is a singleton class, the
    ///constructor is set private to avoid instantiation of the
    ///object.
    Logger();

    ///The map holding the data corresponding the Record() method.
    static map<string, LogEvent> PrflMap;

    ///The stack used by Tic() and Toc().
    static stack<double> TicStack;

    ///The stack to have cumulative 
    static stack<double> FlopStack;
    
    ///The file name for the logger.
    static string log_file;
};

//Profiling Macros
#ifdef PROFILING
#define PROFILESTART() (Logger::Tic())
#define PROFILEEND(str,flps) (                                  \
        Logger::Record(__FUNCTION__, str, Logger::Toc(), flps))
#define PROFILECLEAR() (Logger::PurgeProfileHistory())
#define PROFILEREPORT(format) (Logger::Report(format))
#else
#define PROFILESTART()
#define PROFILECLEAR()
#define PROFILEEND(str,flps)
#define PROFILEREPORT(format)
#endif //PROFILING

//Debugging macros
#ifndef NDEBUG 
#define LOG(msg) (Logger::Log(msg))
#else
#define LOG(msg) 
#endif //NDEBUG

//Printing macro
#define PRINTFRMT scientific<<setprecision(4)
#define CERR(str,endline,action) (                                      \
        std::cerr<<endl<<str                                            \
        <<"\n         File     : "<< __FILE__                           \
        <<"\n         Line     : "<<__LINE__                            \
        <<"\n         Function : "<<__FUNCTION__<<endline,              \
        action                                                          \
                                                                        )

#ifndef NDEBUG
#define COUTDEBUG(str) (std::cout<<str)
#else
#define COUTDEBUG(str)
#endif //NDEBUG
                                                                        
#ifdef VERBOSE
#define COUT(str) (std::cout<<str)
#else
#define COUT(str)
#endif //VERBOSE

//Timing macro
#define GETSECONDS()(omp_get_wtime())

// #define GETSECONDS()(                                       \
//         gettimeofday(&tt, &ttz),                            \
//         (double)tt.tv_sec + (double)tt.tv_usec / 1000000.0)

#endif //_LOGGER_H_
