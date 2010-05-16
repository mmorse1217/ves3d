#include "VectorsTest.h"

int main(int argc, char* argv[])
{
//     VectorsTest<float,CPU> vectest;
//     vectest.PerformAll();

    Device<CPU> dev;
    int p(12);
    int num_vecs(3);
    Vectors<float,CPU> vec(&dev, p, num_vecs);
    Scalars<float,CPU> sc(&dev, p, num_vecs);

    typename Scalars<float,CPU>::iterator it;
    for(it = vec.begin(); it != vec.end(); ++it)
        *it = 2;
    
    DotProduct(vec,vec,sc);
    cout<<sc<<endl;
    
//     cout<<vec<<endl;
    
    return 0;
}

    


