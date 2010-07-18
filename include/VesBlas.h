#ifndef _VESBLAS_H
#define _VESBLAS_H
   
#ifdef HAS_MKL
  #include "HasMkl.h"

#elif  HAS_ATLAS
  #include "HasAtlas.h"

#elif  HAS_BLAS
  #include "HasBlas.h"

#endif

#ifdef GPU_ACTIVE
#include "cublas.h"
#endif //GPU_ACTIVE

#endif //_VESBLAS_H_

