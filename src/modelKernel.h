//
//  BayzR --- modelKernel.h
//
//  Derived class of modelMatrix with small differences:
//   - a 'Kernel' is a similarity matrix stored as eigenvectors and values
//   - the M pointer inherits from Matrix and points to eigenvectors on which the same methods
//     as from matrix can be used. But M is allocated as a dataKernel object so M can hold more.
//     (dataKernel derives from dataMatrix so the same pointer can be used).
//   - additional information in dataKernel is eval, nEvalUsed, ...
//   - correct & collect methods can be used from Matrix, but sample() needs re-implementation
//

#ifndef modelKernel_h
#define modelKernel_h

#include <Rcpp.h>
#include <cmath>
#include "modelMatrix.h"
#include "dataKernel.h"
#include "dataFactor.h"

class modelKernel : public modelMatrix {
   
public:
   
   modelKernel(Rcpp::DataFrame &d, size_t col) : modelMatrix(d, col) {
      Rcpp::RObject col_asRObject = d[col];
      M = new dataKernel(col_asRObject);
      F = new dataFactor(d, col);
   }
   
   ~modelKernel() {
      delete M;
      delete F;
   }
   
   // correct & collect methods not needed here, the ones inherited from modelMatrix are OK

   // need to finish sample ...
   void sample() {
      
   }
   
};

#endif /* modelKernel_h */
