#include<iostream>
#include<exception>
#include "SHVectors.h"

using namespace std;

int main(int argc, char* argv[])
{
    cout<<" ---------------------------"<<endl;
    cout<<" Testing the SHVectors class"<<endl;
    cout<<" ---------------------------"<<endl<<endl;
    

    cout<<" Default constructor and destructor."<<endl;
    cout<<" -----------------------------------"<<endl;
    {
        SHVectors<double> sf;
        sf.p_ = 4;
        sf.number_of_vectors_ = 5;
        sf.number_of_functions_ = 15;
        
        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_vectors    [5]: "<<sf.number_of_vectors_<<endl;
        cout<<" number_of_functions [15]: "<<sf.number_of_functions_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [600]: "<<sf.GetDataLength()<<endl;

	int dLen(1);
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
        SHVectors<float> sf(4,5);

        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_functions [15]: "<<sf.number_of_functions_<<endl;
        cout<<" number_of_vectors    [5]: "<<sf.number_of_vectors_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [600]: "<<sf.GetDataLength()<<endl;
        cout<<" GetVecLength()     [120]: "<<sf.GetVecLength()<<endl;
        
        int dLen(sf.GetDataLength()), idx;
        float *data_in = new float[dLen];
        for(idx=0;idx<dLen;++idx)
            data_in[idx] = idx;
        
        try
        {
            sf.SetData(data_in);
        } 
        catch(range_error err)
        {
            cerr<<" SetData()-data_      [0]: "<<err.what()<<endl;
        }
        delete[] data_in;
        cout<<*(sf.data_+idx-1)-idx+1<<endl<<endl;

    }

    cout<<" Constructor with size and data input."<<endl;
    cout<<" -------------------------------------"<<endl;
    {
        int p(4), nVec(7);
        int dLen(6*p*(p+1)*nVec);
        double *data_in = new double[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        SHVectors<double> sf(p,nVec,data_in);
        delete[] data_in;

        cout<<" p                    [4]: "<<sf.p_<<endl;
        cout<<" number_of_vectors    [7]: "<<sf.number_of_vectors_<<endl;
        cout<<" GetFunLength()      [40]: "<<sf.GetFunLength()<<endl;
        cout<<" GetDataLength()    [840]: "<<sf.GetDataLength()<<endl;
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
            sf.GetFunctionAt(3*nVec);
        }
        catch(range_error err)
        {
            cerr<<" Function index range    : "<<err.what()<<endl;
        }

        const double *fp=sf.GetFunctionAt(1);
        sf.SetFunctionAt(fp,0);
        cout<<" SetFunctionAt()      [0]: "<<*sf.GetFunctionAt(0)-*sf.GetFunctionAt(1)<<endl;
    }
    
    {
        int p(8), nVec(2);
        int fLen(2*p*(p+1));
        int dLen(6*p*(p+1)*nVec);
        double *data_in = new double[dLen];
        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = ii;

        SHVectors<double> sf(p,nVec,data_in),sf2(p,nVec);
        SHScalars<double> dp(p,nVec);
        
        DotProduct(sf,sf,dp);
        cout<<" Dot product       [5 50]: "<<
            *dp.GetFunctionAt(0)<<", "<<*(dp.GetFunctionAt(1)+fLen-1)<<endl;

        CrossProduct(sf,sf,sf2);
        cout<<" Cross product      [0 0]: "<<
            *sf2.GetFunctionAt(0)<<", "<<*(sf2.GetFunctionAt(5)+fLen-1)<<endl;
        
        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = -ii;
        sf2.SetData(data_in);
        CrossProduct(sf,sf2,sf2);
        cout<<" Cross product      [0 0]: "<<
            *sf2.GetFunctionAt(0)<<", "<<*(sf2.GetFunctionAt(5)+fLen-1)<<endl;
    
        for(int ii=0;ii<dLen;++ii)
            data_in[ii] = 1;
        sf.SetData(data_in);

        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = ii+3;
        sf2.SetData(data_in);

        CrossProduct(sf,sf2,sf2);
        cout<<" Cross product   [1 -2 1]: "<<
            *sf2.GetFunctionAt(0)<<", "<<*sf2.GetFunctionAt(1)<<", "<<*sf2.GetFunctionAt(2)<<endl;
        delete[] data_in;
    }

    {
        int p(8), nFun(2), nVec(2);
        int fLen(2*p*(p+1));
        int dLen(2*p*(p+1)*nFun);
        
        double *data_in = new double[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        SHScalars<double> sf(p,nFun,data_in);
        delete[] data_in;

        dLen =6*p*(p+1)*nVec;
        data_in = new double[dLen];
        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = 1;

        SHVectors<double> sv(p,nVec,data_in),sv2(p,nVec);


        AxPy(2.0,sv,0.0,sv2);
        cout<<" (inherited) AxPy()     [2]: "<<sv2.data_[0]<<endl;
        
        AxPy(sf,sv,0.0,sv);
        cout<<" AxPy()                 [5]: "<<sv.data_[5]<<endl;

        xDy(sv,sf,sv);
        cout<<" AxPy()                 [1]: "<<sv.data_[5]<<endl;
        
        AxPy(sf,sv,sv2,sv2);
        cout<<" AxPy()                 [7]: "<<sv2.data_[5]<<endl;

    }

    cout<<" Memory allocation."<<endl;
    cout<<" -----------------"<<endl;
    {
        int p(64), nVec(200);
        int dLen(6*p*(p+1)*nVec);
        double *data_in = new double[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        for(int idx=0;idx<500;++idx)
        {
            SHVectors<double> sf(p,nVec,data_in);
        }

        delete[] data_in;
    }

    
    cout<<" ------------ "<<endl;
    cout<<" End of test. "<<endl;
    cout<<" ------------ "<<endl;
}



