#include "Logger.h"

stack<double> Logger::TicStack;
map<string, LogEvent> Logger::PrflMap;
string Logger::log_file;

//The variables for the timing macro.
// struct timeval  tt;
// struct timezone ttz;

template<typename T>
void PrintLogEvent(const pair<T, LogEvent> &ev)
{
    int print_length = 27;
    string printstr = ev.second.fun_name;
    printstr.resize(print_length,' ');
    printstr= "     " + printstr +
        "%-10u \t %-4.3e \t  %-4.3e  \t %-4.3e\n";

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
    return(now);
}

double Logger::Toc()
{
    if(Logger::TicStack.empty())
        CERR(" There is no matching Logger::Tic() call.",endl,sleep(0));

    double toc = Logger::Now();
    toc -= Logger::TicStack.top();
    Logger::TicStack.pop();

    return(toc);
}

void Logger::Record(string fun_name, string prefix, double time, double flop)
{
    fun_name = prefix+fun_name;
#pragma omp critical
    if(Logger::PrflMap.count(fun_name))
    {
        Logger::PrflMap[fun_name].num_calls++;
        Logger::PrflMap[fun_name].time +=time;
        Logger::PrflMap[fun_name].flop +=flop;
    }
    else
    {
        LogEvent ev;
        ev.fun_name = fun_name;
        ev.num_calls = 1;
        ev.time = time;
        ev.flop = flop;
        Logger::PrflMap.insert(make_pair(fun_name, ev));
    }
}

void Logger::Report(enum ReportFormat rf) 
{

    multimap<double, LogEvent> ReportMap;
    map<string, LogEvent>::iterator it;

    for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
        it->second.flop_rate = it->second.flop/it->second.time/1e9;

    switch ( rf )
    {
        case SortFunName:
            COUT("\n =========================================================================================="<<endl
                <<"   >Function name               Calls            Total time         GFlop         GFlop/sec"
                <<"\n ------------------------------------------------------------------------------------------"<<endl);
            for_each(PrflMap.begin(), PrflMap.end(), &PrintLogEvent<string>);
            break;

        case SortNumCalls:
            COUT("\n ==========================================================================================="<<endl
                <<"    Function name              >Calls            Total time         GFlop        GFlop/sec"
                <<"\n -------------------------------------------------------------------------------------------"<<endl);
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(make_pair(it->second.num_calls, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

        case SortTime:
            COUT("\n ==========================================================================================="<<endl
                <<"    Function name               Calls           >Total time         GFlop        GFlop/sec"
                <<"\n -------------------------------------------------------------------------------------------"<<endl);
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(make_pair(it->second.time, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;
            
        case SortFlop:
            COUT("\n ==========================================================================================="<<endl
                <<"    Function name               Calls            Total time        >GFlop        GFlop/sec"
                <<"\n -------------------------------------------------------------------------------------------"<<endl);
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(make_pair(it->second.flop, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;

        case SortFlopRate:
            COUT("\n ==========================================================================================="<<endl
                <<"    Function name               Calls            Total time         GFlop       >GFlop/sec"
                <<"\n -------------------------------------------------------------------------------------------"<<endl);
            for (it = Logger::PrflMap.begin();it != Logger::PrflMap.end(); ++it)
                ReportMap.insert(make_pair(it->second.flop_rate, it->second));
            for_each(ReportMap.begin(), ReportMap.end(), &PrintLogEvent<double>);
            break;
            
    }
    COUT(" ==========================================================================================="<<endl);
    
    if(TicStack.size())
        CERR("\n There may be unbalanced Tic() and Toc() calls.",endl,sleep(0));
}

void Logger::SetLogFile(string file_name)
{
    log_file = file_name;
}
    
void Logger::Log(const char *event)
{ 
    ofstream fl(log_file.c_str(), ios::app);

    if(!fl)
    {
        CERR(" Could not open the log file." <<endl
            <<" File name : " << log_file <<"."<<endl
            <<"\n Log event: "<<event,endl,exit(1));
    }
    else
    {
        fl<<event<<endl;
    }
    
    fl.close();
}
