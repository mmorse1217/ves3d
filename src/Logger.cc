#include "Logger.h"

#include <algorithm> //to get for_each()
#include <fstream>   //ofstream type
#include <unistd.h>  //to get sleep()


std::stack<double> Logger::TicStack;
std::stack<double> Logger::FlopStack;
std::map<std::string, LogEvent> Logger::PrflMap;
std::string Logger::log_file;

// //The variables for the timing macro.
// struct timeval  tt;
// struct timezone ttz;

template<typename T>
void PrintLogEvent(const std::pair<T, LogEvent> &ev)
{
    int print_length = 27;
    std::string printstr = ev.second.fun_name;
    printstr.resize(print_length,' ');
    printstr= "  " + printstr +
        "%-8u \t %-4.3e \t %-4.3e \t %-4.3e \n";

#ifdef VERBOSE
    printf(printstr.c_str(), ev.second.num_calls, ev.second.time,
        ev.second.flop/1e9, ev.second.flop_rate);
#endif
}

double Logger::Now()
{
    return(GETSECONDS());
}

double Logger::Tic()
{
    double now = Logger::Now();
    Logger::TicStack.push(now);

    if ( Logger::FlopStack.empty() )
        Logger::FlopStack.push(0); //in order to have the total number of flops

    Logger::FlopStack.push(0);
    return(now);
}

double Logger::Toc()
{
    double toc;
    if(Logger::TicStack.empty())
        CERR("There is no matching Logger::Tic() call.","", toc=0);
    else
    {
        toc = Logger::Now();
        toc -= Logger::TicStack.top();
        Logger::TicStack.pop();
    }
    return(toc);
}

void Logger::Record(std::string fun_name, std::string prefix,
    double time, double flop)
{
#pragma omp critical (loggerRecord)
    {
        fun_name = prefix+fun_name;
        flop += Logger::FlopStack.top();
        Logger::FlopStack.pop();
        if ( !Logger::FlopStack.empty() )
            Logger::FlopStack.top() += flop;

        if(Logger::PrflMap.count(fun_name))
        {
            Logger::PrflMap[fun_name].num_calls++;
            Logger::PrflMap[fun_name].time +=time;
            Logger::PrflMap[fun_name].flop +=flop;
        }
        else
        {
            LogEvent ev;
            ev.fun_name  = fun_name;
            ev.num_calls = 1;
            ev.time      = time;
            ev.flop      = flop;
            Logger::PrflMap.insert(make_pair(fun_name, ev));
        }
    }
}

void Logger::PurgeProfileHistory()
{
    Logger::PrflMap.clear();

    while ( !Logger::TicStack.empty() )
        Logger::TicStack.pop();

    while ( !Logger::FlopStack.empty() )
        Logger::FlopStack.pop();
}

void Logger::Report(enum ReportFormat rf)
{
    std::multimap<double, LogEvent> ReportMap;
    std::map<std::string, LogEvent>::iterator it;

    for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
        it->second.flop_rate = it->second.flop/it->second.time/1e9;

    switch ( rf )
    {
        case SortFunName:
            COUT("=====================================================================================\n"
                <<" >Function name              Calls       Total time      GFlop           GFlop/sec\n"
                <<"------------------------------------------------------------------------------------");
            for_each(PrflMap.begin(), PrflMap.end(), &PrintLogEvent<std::string>);
            break;

        case SortNumCalls:
            COUT("=====================================================================================\n"
                <<"  Function name             >Calls       Total time      GFlop           GFlop/sec\n"
                <<"------------------------------------------------------------------------------------");
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(std::make_pair(it->second.num_calls, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

        case SortTime:
            COUT("=====================================================================================\n"
                <<"  Function name              Calls      >Total time      GFlop           GFlop/sec\n"
                <<"------------------------------------------------------------------------------------");
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(std::make_pair(it->second.time, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

        case SortFlop:
            COUT("=====================================================================================\n"
                <<"  Function name              Calls       Total time     >GFlop           GFlop/sec\n"
                <<"------------------------------------------------------------------------------------");
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(std::make_pair(it->second.flop, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

        case SortFlopRate:
            COUT("=====================================================================================\n"
                <<"  Function name              Calls       Total time      GFlop          >GFlop/sec\n"
                <<"------------------------------------------------------------------------------------");
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(std::make_pair(it->second.flop_rate, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

    }
    COUT("=====================================================================================");

    if(TicStack.size())
        CERR("There may be unbalanced Tic() and Toc() calls.","",sleep(0));
}

void Logger::SetLogFile(std::string file_name)
{
    log_file = file_name;
}

void Logger::Log(const char *event)
{
    std::ofstream fl(log_file.c_str(), std::ios::app);

    if(!fl)
    {
        CERR("Could not open the log file." <<std::endl
            <<" File name : " << log_file <<"."<<std::endl
            <<"\n Log event: "<<event,"",exit(1));
    }
    else
    {
        fl<<event<<std::endl;
    }

    fl.close();
}

double Logger::GetFlops()
{
    if ( !Logger::FlopStack.empty() )
        return(Logger::FlopStack.top());
    else
        return(0);
}

std::ostream& alert(std::ostream& os)
{
    if(os.iword(alert_xalloc) == false)
        os<<"\033[1;31m";
    else
        os<<"\033[0m";

    os.iword(alert_xalloc) = !os.iword(alert_xalloc);
    return os;
}

std::ostream& emph(std::ostream& os)
{
    if(os.iword(emph_xalloc) == false)
        os<<"\033[1;36m";
    else
        os<<"\033[0m";

    os.iword(emph_xalloc) = !os.iword(emph_xalloc);
    return os;
}
