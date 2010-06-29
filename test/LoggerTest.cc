#include "Logger.h"
#include <fstream>

using namespace std;

void Wait(int sec)
{
    PROFILESTART();
    LOG("iWait");
    sleep(sec);
    PROFILEEND("i",0);
}

void Wait(double sec)
{
    PROFILESTART();
    LOG("dWait");
    sleep(sec);
    PROFILEEND("d",0);
}

int Dummy(int ii)
{
    PROFILESTART();
    LOG("Dummy");
    ii = ii + 10;
    Wait((double) 0.5);
    PROFILEEND("",0);
    return(ii);
}

int main(int argc, char** argv)
{
    PROFILESTART();
    string log_file("LoggerTest.out");
    Logger::SetLogFile(log_file);

    Logger::Tic();
    int imax(3);
    for(int ii=0;ii<imax;++ii)
        Dummy(ii);
    cout<<" "<<imax<<" calls to Dummy : "<<Logger::Toc()<<endl;
   
    Logger::Tic();
    Wait(1);
    Logger::Tic();
    sleep(1);
    cout<<" Inner Tic/Toc : "<<Logger::Toc()<<endl;
    cout<<" Outer Tic/Toc : "<<Logger::Toc()<<endl;
        
    Wait((double) 1.0);
    PROFILEEND("",0);
    PROFILEREPORT(SortTime);
    
    //printing the log file
    COUT("\n\n  The content of the log file: \n ==============================\n");
    
    ifstream file(log_file.c_str());
    string line;
    while ( getline ( file, line ) )
        cout<<"  "<<line<<endl;

    file.close();

    remove(log_file.c_str());
}
