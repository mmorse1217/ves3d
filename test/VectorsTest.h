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
    Vectors<T> &vec;

    ScalarsTest(Vectors<T> &vecIn) :
        vec(veccIn)
    { }

    bool performAll()
    {
        int p(4), nVec(7);
        int dLen(6*p*(p+1)*nVec);
        T *data_in = new T[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        vec.Resize(p,nVec);
        vec.SetData(data_in);
        delete[] data_in;

        cout<<" GetFunLength()      [40]: "<<vec.GetFunLength()<<endl;
        cout<<" GetDataLength()    [840]: "<<vec.GetDataLength()<<endl;
        cout<<" GetFunctionAt() [0, 200]: "<<*vec.GetFunctionAt(0)<<", "<<*vec.GetFunctionAt(5)<<endl;
 

        const T *fp=vec.GetFunctionAt(1);
        vec.SetFunctionAt(fp,0);
        cout<<" SetFunctionAt()      [0]: "<<*vec.GetFunctionAt(0)-*vec.GetFunctionAt(1)<<endl;

        return false;
    }
};
        Vectors<T> vec(p,nVec,data_in),vec2(p,nVec);
        Scalars<T> dp(p,nVec);
        
        DotProduct(vec,vec,dp);
        cout<<" Dot product       [5 50]: "<<
            *dp.GetFunctionAt(0)<<", "<<*(dp.GetFunctionAt(1)+fLen-1)<<endl;

        CrossProduct(vec,vec,vec2);
        cout<<" Cross product      [0 0]: "<<
            *vec2.GetFunctionAt(0)<<", "<<*(vec2.GetFunctionAt(5)+fLen-1)<<endl;
        
        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = -ii;
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
    }

    {
        int p(8), nFun(2), nVec(2);
        int fLen(2*p*(p+1));
        int dLen(2*p*(p+1)*nFun);
        
        T *data_in = new T[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;

        Scalars<T> vec(p,nFun,data_in);
        delete[] data_in;

        dLen =6*p*(p+1)*nVec;
        data_in = new T[dLen];
        for(int ii=0;ii<3*nVec;++ii)
            for(int idx=0;idx<fLen;++idx)
                data_in[ii*fLen+idx] = 1;

        Vectors<T> sv(p,nVec,data_in),sv2(p,nVec);


        AxPy(2.0,sv,0.0,sv2);
        cout<<" (inherited) AxPy()     [2]: "<<sv2.data_[0]<<endl;
        
        AxPy(vec,sv,0.0,sv);
        cout<<" AxPy()                 [5]: "<<sv.data_[5]<<endl;

        xDy(sv,vec,sv);
        cout<<" AxPy()                 [1]: "<<sv.data_[5]<<endl;
        
        AxPy(vec,sv,sv2,sv2);
        cout<<" AxPy()                 [7]: "<<sv2.data_[5]<<endl;

    }

//     cout<<" Memory allocation."<<endl;
//     cout<<" -----------------"<<endl;
//     {
//         int p(64), nVec(200);
//         int dLen(6*p*(p+1)*nVec);
//         T *data_in = new T[dLen];
//         for(int idx=0;idx<dLen;++idx)
//             data_in[idx] = idx;

//         for(int idx=0;idx<500;++idx)
//         {
//             Vectors<T> vec(p,nVec,data_in);
//         }

//         delete[] data_in;
//     }

    
    cout<<" ------------ "<<endl;
    cout<<" End of test. "<<endl;
    cout<<" ------------ "<<endl;
}



