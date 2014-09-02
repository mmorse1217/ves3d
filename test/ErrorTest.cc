#include "Error.h"
#include "Logger.h"

#include <sstream>
#include <cassert>

using namespace std;

#ifndef Doxygen_skip

Error_t local_cb_a(const Error_t &err)
{
    COUT("CB_A: Received "<<err<<endl);
    return err;
}

Error_t local_cb_b(const Error_t &err)
{
    COUT("CB_B: Received "<<err<<endl);
    return err;
}
#endif //Doxygen_skip

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Error Test:"
        <<"\n ==============================\n");

    CHK(ErrorEvent::Success);
    CHK(ErrorEvent::InvalidParameter);

    SET_ERR_CALLBACK(&local_cb_a);
    CHK(ErrorEvent::InvalidDevice);
    CHK_CB(ErrorEvent::SetOnActiveDevice,&local_cb_b);

    if (ERRORSTATUS())
        COUT("There was some error"<<endl);

    PRINTERRORLOG();
    CLEARERRORHIST();

    int i;
    for(i=0;i<10;++i){
        if(i==3)
            CHK(ErrorEvent::SolverDiverged);
	COUT(i<<endl);
        BREAKONERROR();
    }
    COUT("Stopped at i="<<i<<endl);

    {
        ostringstream stream;
        string str;
        stream<<ErrorEvent::Success;
        COUT("Overloaded streaming operator"<<stream.str()<<endl);
        assert(stream.str()=="Success");
    }

    PRINTERRORLOG();
}
