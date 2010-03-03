/**
 * @file   ScalarsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include "Scalars.h"
#include <iostream>

using namespace std;

template<typename T>
class ScalarsTest
{
  public:
    Scalars<T> &sc;

    ScalarsTest(Scalars<T> &scIn) :
        sc(scIn)
    { }

    bool performAll()
    {
        cout<<" ---------------------------"<<endl;
        cout<<" Testing the SHScalars class"<<endl;
        cout<<" ---------------------------"<<endl<<endl;

        
        int p(4), nFun(7);
        int dLen(2*p*(p+1)*nFun);
        T *data_in = new T[dLen];
        for(int idx=0;idx<dLen;++idx)
            data_in[idx] = idx;
        
        sc.Resize(p,nFun);
        sc.SetData(data_in);
            
        cout<<" GetFunLength()      [40]: "<<sc.GetFunLength()<<endl;
        cout<<" GetDataLength()    [280]: "<<sc.GetDataLength()<<endl;
        cout<<" GetFunctionAt() [0, 200]: "<<*sc.GetFunctionAt(0)<<", "<<*sc.GetFunctionAt(5)<<endl;
            
        const T *fp=sc.GetFunctionAt(1);
        sc.SetFunctionAt(fp,0);
        cout<<" SetFunctionAt()       [0]: "<<*sc.GetFunctionAt(0)-*sc.GetFunctionAt(1)<<endl;
        
        axpb((T) 1.0,sc,(T) 2.0,sc);
        cout<<" axpb()               [42]: "<<*sc.GetFunctionAt(0)<<endl;
        axpy((T) 1.0,sc,sc,sc);
        cout<<" axpy()               [84]: "<<*sc.GetFunctionAt(0)<<endl;
        axpb((T) .5,sc,(T) -2.0,sc);
        cout<<" axpb()               [40]: "<<*sc.GetFunctionAt(0)<<endl;
        xy(sc,sc,sc);
        cout<<" xy()               [1600]: "<<*sc.GetFunctionAt(0)<<endl;
        sc.Sqrt();
        cout<<" Sqrt()            [40,41]: "<<*sc.GetFunctionAt(0)<<", "<<*(sc.GetFunctionAt(0)+1)<<endl<<endl;
        xyInv(sc,sc,sc);
        cout<<" xyInv()               [1]: "<<*sc.GetFunctionAt(0)<<endl<<endl;
            
            
        cout<<" ------------ "<<endl;
        cout<<" End of test. "<<endl;
        cout<<" ------------ "<<endl;
        
        return false;
    }
};

    


