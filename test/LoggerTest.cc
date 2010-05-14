#include "Logger.h"

using namespace std;

int Dummy(int ii)
{
    PROFILESTART();
    LOG("Dummy");
    ii = ii + 10;
    PROFILEEND("0",0);
    return(ii);
}

void Wait(int sec)
{
    PROFILESTART();
    LOG("Wait");
    sleep(sec);
    PROFILEEND("1",0);
}

void Wait(double sec)
{
    PROFILESTART();
    LOG("Wait");
    sleep(sec);
    PROFILEEND("2",0);
}

int main(int argc, char** argv)
{
    PROFILESTART();
    //Logger::SetLogFile("log");

    Logger::Tic();
    for(int jj=0;jj<1e3;++jj)
        for(int ii=0;ii<1e3;++ii)
            Dummy(ii);
    cout<<Logger::Toc()<<endl;
   
    Logger::Tic();
    Wait(1);
    Logger::Tic();
    sleep(1);
    cout<<Logger::Toc()<<endl;
    cout<<Logger::Toc()<<endl;
    
    Wait((double) 1.0);
    PROFILEEND("",0);
    PROFILEREPORT(SortTime);
}
