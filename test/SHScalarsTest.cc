#include<iostream>
#include<exception>
#include "SHScalars.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout<<" ---------------------------"<<endl;
    cout<<" Testing the SHScalars class"<<endl;
    cout<<" ---------------------------"<<endl<<endl;
    

    cout<<" Default constructor and destructor."<<endl;
    cout<<" -----------------------------------"<<endl;
    {
        SHScalars<double> sf;
        sf.p_ = 4;
        sf.number_of_functions_ = 5;
        
        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_functions  [5]: "<<sf.number_of_functions_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [200]: "<<sf.GetDataLength()<<endl;
        
	int dLen(200);
        double *data_in = new double[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;
        try{
            sf.SetData(data_in);
        } catch(range_error err) {
            cerr<<" data_ array allocation  :"<<err.what()<<endl;
        }
        delete[] data_in;
        cout<<endl<<endl;
    }

    cout<<" Constructor with size input."<<endl;
    cout<<" ----------------------------"<<endl;
    {
        SHScalars<float> sf(4,5);

        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_functions  [5]: "<<sf.number_of_functions_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [200]: "<<sf.GetDataLength()<<endl;
        cout<<" SetData()-data_      [0]: ";

        int dLen(200), idx;
        float *data_in = new float[dLen];
        for(idx=0;idx<dLen;++idx)
            data_in[idx] = idx;
        
        try
        {
            sf.SetData(data_in);
        } 
        catch(range_error err)
        {
            cerr<<err.what()<<endl;
        }
        delete[] data_in;
	cout<<*(sf.data_+idx-1)-idx+1<<endl<<endl;

    }

    cout<<" Constructor with size and data input."<<endl;
    cout<<" -------------------------------------"<<endl;
    {
        int p(4), nFun(7);
        int dLen(2*p*(p+1)*nFun);
        double *data_in = new double[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        SHScalars<double> sf(p,nFun,data_in);
        delete[] data_in;

        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_functions  [7]: "<<sf.number_of_functions_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [280]: "<<sf.GetDataLength()<<endl;
        cout<<" GetFunctionAt() [0, 200]: "<<*sf.GetFunctionAt(0)<<", "<<*sf.GetFunctionAt(5)<<endl;
 
        try
        {
            sf.GetFunctionAt(-1);
        }
        catch(range_error err)
        {
	  cerr<<" Function index range    : "<<err.what()<<endl;
        }

        
        try
        {
	  sf.GetFunctionAt(7);
        }
        catch(range_error err)
        {
	  cerr<<" Function index range    : "<<err.what()<<endl;
        }

        const double *fp=sf.GetFunctionAt(1);
        sf.SetFunctionAt(fp,0);
        cout<<" SetFunctionAt()       [0]: "<<*sf.GetFunctionAt(0)-*sf.GetFunctionAt(1)<<endl<<endl;

        cout<<" ------------ "<<endl;
        cout<<" End of test. "<<endl;
        cout<<" ------------ "<<endl;
        
    }
}


