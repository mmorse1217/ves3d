#include "Logger.h"
#include <fstream>

using namespace std;

#ifndef Doxygen_skip

void Wait(int sec)
{
    PROFILESTART();
    LOG("Wait (integer arg)");
    sleep(sec);
    PROFILEEND("i",0);
}

void Wait(double sec)
{
    PROFILESTART();
    LOG("Wait (double arg)");
    sleep((unsigned int) sec);
    PROFILEEND("d",1000);
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

#endif //Doxygen_skip

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Logger Test:"
        <<"\n ==============================\n");
    sleep(1);

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
    COUT("\n\n  Content of the log file: \n"
        <<" ==============================\n");

    ifstream file(log_file.c_str());
    string line;
    while ( getline ( file, line ) )
        cout<<"  "<<line<<endl;

    file.close();

    remove(log_file.c_str());
    sleep(1);
}
