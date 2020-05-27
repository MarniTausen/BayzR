//
//  modelMatrix.h
//  rbayz
//
//  Defines the computational methods when design matrix is a (real) matrix:
//    -> has pointer to dataMatrix object
//    -> has pointer to dataFactor object
//  and defines residual de/correct and collect lhs/rhs methods working on this
//  kind of objects.
//  This is not yet a concrete class, derived classes differ mostly in the
//  constructors that define how different kinds of matrix data is prepared.
//
//  Created by Luc Janss on 30/08/2019.
//

#ifndef modelMatrix_h
#define modelMatrix_h

#include <Rcpp.h>
#include <cmath>
#include "dataMatrix.h"

class modelMatrix : public modelBase {
   
public:
   
   modelMatrix(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      *M = NULL;
      *F = NULL;
      lhs = 0.0l;
      rhs = 0.0l;
   }
   
   ~modelMatrix() {
   }
   
protected:

   void resid_correct(size_t col);
   void resid_decorrect(size_t col);
   void collect_lhs_rhs(size_t col);

   dataMatrix *M;
   dataFactor *F;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   
};

#endif /* modelMatrix_h */
