#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <iostream>
#include <sys/time.h>
#define GetSeconds()   (gettimeofday(&tt, &ttz), \
    (double)tt.tv_sec + (double)tt.tv_usec / 1000000.0)

struct timeval  tt;
struct timezone ttz;

enum LogLevel {NoLogging = 100, 
               FlopCount = 101, 
               Profiling = 102,
               Debug     = 103};

class Logger
{
  public:
    static const enum LogLevel the_log_level_;
    static void Log(enum LogLevel events_log_level, const char *msg = 0);
    static void Add2Flops(double flops_in);

    static double Now();
    static double Tic();
    static double Toc();

    static void TearDown();
  
    static double GetGFlops();
  private:
    Logger();
    static double flop_count_;
    static double last_tic_;
};

void Logger::Log(enum LogLevel events_log_level, const char *msg)
{

    if(events_log_level>the_log_level_)
        return;

//     switch (events_log_level)
//     {
//         case Debug:

//         case Profiling:
            
//         case FlopCount:
//     }
}

void Logger::Add2Flops(double flops_in)
{
    flop_count_ += flops_in;
}

double Logger::GetGFlops()
{
    return(flop_count_/1e9);
}

Logger::Logger() {}

double Logger::Now()
{
    return(GetSeconds());
}

double Logger::Tic()
{
    last_tic_ = Logger::Now();
    return(last_tic_);
}

double Logger::Toc()
{
    double toc = Logger::Now();
    toc -= last_tic_;
    last_tic_ = 0;
    return(toc);
}

void Logger::TearDown() 
{
    if(the_log_level_ != NoLogging)
    {
        cout<<"\n ++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
        cout<<"  The Total number of Operations : "
            <<flop_count_ / 1e9
            <<" GFlops."<<endl;
        cout<<" ++++++++++++++++++++++++++++++++++++++++++++++++"<<endl<<endl;
    }
}


double Logger::flop_count_ = 0;
double Logger::last_tic_ = 0;


#ifndef NDEBUG
const enum LogLevel Logger::the_log_level_ = Debug;
#else

#ifdef NOLOGGING
const enum LogLevel Logger::the_log_level_ = NoLogging;
#elif PROFILING
const enum LogLevel Logger::the_log_level_ = Profiling;
#else 
const enum LogLevel Logger::the_log_level_ = FlopCount;
#endif //Logging level

#endif //NDEBUG

#endif //_LOGGER_H_
