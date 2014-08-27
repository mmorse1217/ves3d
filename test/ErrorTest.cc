#include "Error.h"
#include "Logger.h"

using namespace std;

#ifndef Doxygen_skip

Error_t local_cb_a(const Error_t &err)
{
    COUT("CB_A: Received "<<err<<endl);
}

Error_t local_cb_b(const Error_t &err)
{
    COUT("CB_B: Received "<<err<<endl);
}
#endif //Doxygen_skip

int main(int argc, char** argv)
{
    COUT("\n ==============================\n"
        <<"  Error Test:"
        <<"\n ==============================\n");

    CHK(Success);
    CHK(InvalidParameter);

    SET_ERR_CALLBACK(&local_cb_a);
    CHK(InvalidDevice);
    CHK_CB(SetOnActiveDevice,&local_cb_b);

    if (ERRORSTATUS())
        COUT("There was some error"<<endl);

    PRINTERRORLOG();
    CLEARERRORHIST();

    int i;
    for(i=0;i<10;++i){
        if(i=3)
            CHK(SolverDiverged);

        BREAKONERROR();
    }
    COUT("Stopped at i="<<i<<endl);

    PRINTERRORLOG();
}
