/**
 * @file   VectorsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 20:11:43 2010
 * 
 * @brief  The tester for the Vectors class.
 */


#include "Vectors.h"
#include <iostream>

using namespace std;

template<typename T>
class VectorsTest
{
  public:
    Device<T> &device;

    VectorsTest(Device<T> &device_in) :
        device(device_in)
    { }

    bool performAll()
    {
        int p(4), nVec(2);
        Vectors<T> vec(device,p,nVec);
        int fLen(2*p*(p+1));
        int dLen(6*p*(p+1)*nVec);
        T *data_in = new T[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;
        
        vec.SetData(data_in);
        delete[] data_in;

        cout<<" GetFunLength()      [40]: "<<vec.GetFunLength()<<endl;
        cout<<" GetDataLength()    [840]: "<<vec.GetDataLength()<<endl;
        cout<<" GetFunctionAt() [0, 200]: "<<*vec.GetFunctionAt(0)<<", "<<*vec.GetFunctionAt(5)<<endl;
 

        Vectors<T> vec2(device,p,nVec);
        Scalars<T> dp(device,p,nVec);
        
        DotProduct(vec,vec,dp);
        cout<<" Dot product [8000 80000]: "<<
            *dp.GetFunctionAt(0)<<", "<<*(dp.GetFunctionAt(1))<<endl;

        CrossProduct(vec,vec,vec2);
        cout<<" Cross product      [0 0]: "<<
            *vec2.GetFunctionAt(0)<<", "<<*(vec2.GetFunctionAt(5)+fLen-1)<<endl;
        
        data_in = new T[dLen];
        for(int ii=0;ii<vec.GetDataLength();++ii)
                data_in[ii] = -vec.data_[ii];
        vec2.SetData(data_in);
        CrossProduct(vec,vec2,vec2);
        cout<<" Cross product      [0 0]: "<<
            *vec2.GetFunctionAt(0)<<", "<<*(vec2.GetFunctionAt(5)+fLen-1)<<endl;
    
        for(int ii=0;ii<dLen;++ii)
            data_in[ii] = 1;
        vec.SetData(data_in);

        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = ii+3;
        vec2.SetData(data_in);

        CrossProduct(vec,vec2,vec2);
        cout<<" Cross product   [1 -2 1]: "<<
            *vec2.GetFunctionAt(0)<<", "<<*vec2.GetFunctionAt(1)<<", "<<*vec2.GetFunctionAt(2)<<endl;
        delete[] data_in;

        axpb((T) 2.0,vec, (T) 0.0,vec2);
        cout<<" (inherited) AxPy()   [2]: "<<vec2.data_[0]<<endl;

        dp.SetData(vec2.data_);
        
        xvpw(dp,vec,vec2,vec);
        cout<<" xvpw()               [4]: "<<vec.data_[1]<<endl;

        xvpb(dp,vec, (T) 3,vec);
        cout<<" xvpw()              [11]: "<<vec.data_[1]<<endl;

        return false;
    }
};


