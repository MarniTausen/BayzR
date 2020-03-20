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

class modelMatrix : public modelTerm {
   
public:
   
   dataMatrix(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      *M = NULL;
      *F = NULL;    // make this one derive from modelFactor? - then modelFactor
      lhs = 0.0l;   // constructor can handle factor set-up.
      rhs = 0.0l;
   }
   
   ~dataMatrix() {
   }
   
protected:

   void resid_correct(size_t col) {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] -= par[col] * M->data(F->data(obs),col);
   }
   
   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] += par[col] * M->data(F->data(obs),col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t rowlevel;
      for (size_t obs=0; obs < F->data.size(); obs++) {
         rowlevel = F->data(obs);
         rhs += M->data(rowlevel,col) * residPrec[obs] * resid[obs];
         lhs += M->data(rowlevel,col) * M->data(rowlevel,col) * residPrec[obs];
      }
   }

   dataMatrix *M;
   dataFactor *F;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   
};

#endif /* modelTerm_realmat_h */
