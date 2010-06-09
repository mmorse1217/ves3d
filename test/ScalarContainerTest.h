/**
 * @file   ScalarsTest.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Mar  2 09:55:22 2010
 *
 * @brief  Tester class for the Scalars class.
 */

#include "Logger.h"
#include <iostream>

using namespace std;

template<typename Container>
class ScalarsTest
{
  public:
    bool PerformAll();
};

template<typename Container>
bool ScalarsTest<Container>::PerformAll()
{
    int p = 12;
    int num_funs(1000);

//     Container sc;
//     sc.resize(num_funs, p);
    
//     typename Container::iterator it;
//     for(it = sc.begin(); it != sc.end(); ++it)
//         *it = 10;

//     cout<<sc<<endl;
//     sc.resize(0); cout<<sc<<endl;
//     sc.resize(3); cout<<sc<<endl;

//     Container sc2;
//     sc2.replicate(sc);
//     cout<<sc2<<endl;

//     cout<<sc2.getTheDim()<<" "<<sc2.getNumSubs()<<endl;
    
    int N = 100000;
    Logger::Tic();  
    for(int ii=0; ii<N; ++ii)
    {
        Container sc(num_funs, p);
        sc[0] = sc[2];
}
    cout<<"Creation-destruction time : "<<Logger::Toc()<<endl;

    Container sc(num_funs, p);
    Logger::Tic();  
    for(int ii=0; ii<N; ++ii)
    {
        sc.resize(0);
        sc.resize(num_funs);
        sc[0] = sc[2];
    }
    cout<<"Freeing-allocation time : "<<Logger::Toc()<<endl;

    Logger::Tic();  
    for(int ii=0; ii<N; ++ii)
    {
        Container* sc = new Container(num_funs, p);
        delete sc;
    }
    cout<<"Creation-destruction (new) time : "<<Logger::Toc()<<endl;


    return true;
}
    


