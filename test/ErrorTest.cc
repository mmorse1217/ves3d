#include "Error.h"
#include "Logger.h"
#include "TestTools.h"

#include <sstream>
#include <cassert>

using namespace std;

#ifndef Doxygen_skip

Error_t local_cb_a(const ErrorEvent &err)
{
    COUT("CB_A: Received "<<err);
    return err.err_;
}

Error_t local_cb_b(const ErrorEvent &err)
{
    COUT("CB_B: Received "<<err);
    return err.err_;
}
#endif //Doxygen_skip

int main(int argc, char** argv)
{
    VES3D_INITIALIZE(&argc,&argv,NULL,NULL);
    COUT("\n ==============================\n"
        <<"  Error Test:"
        <<"\n ==============================\n");

    CHK(ErrorEvent::Success);
    CHK(ErrorEvent::InvalidParameterError);

    SET_ERR_CALLBACK(&local_cb_a);
    CHK(ErrorEvent::InvalidDeviceError);
    CHK_CB(ErrorEvent::SetOnActiveDeviceError,&local_cb_b);

    if (ERRORSTATUS())
        COUT("There was some error");

    PRINTERRORLOG();
    CLEARERRORHIST();

    int i;
    for(i=0;i<10;++i){
        if(i==3)
            CHK(ErrorEvent::DivergenceError);
	COUT(i);
        BREAKONERROR();
    }
    COUT("Stopped at i="<<i);

    {
        ostringstream stream;
        string str;
        stream<<ErrorEvent::Success;
        COUT("Overloaded streaming operator"<<stream.str());
	testtools::AssertTrue(stream.str()=="Success","correct string", "bad string");
    }

    PRINTERRORLOG();
    COUT(emph<<"** ErrorTest passed **"<<emph<<std::endl);
    VES3D_FINALIZE();
}
