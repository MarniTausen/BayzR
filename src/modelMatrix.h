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
#include "dataFactor.h"
#include "dataMatrix.h"
#include "modelBase.h"

class modelMatrix : public modelBase {
   
public:
   
   modelMatrix(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      Rcpp::RObject col_asRObject = d[col];
      M = new dataMatrix(col_asRObject);
      F = new dataFactor(d, col);
      lhs = 0.0l;
      rhs = 0.0l;
   }
   
   ~modelMatrix() {
      delete M;
      delete F;
   }
   
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

   void sample() {
      for(size_t k=0; k < M->nColUsed; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         lhs = lhs + (1.0/( M->weights[k]*hpar[0] ));   // lhs with variance added
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));  // Note weights still stored as variances, not inverse
         resid_correct(k);
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->nColUsed; k++)
         ssq += par[k]*par[k]/M->Weights[k];
      hpar[0] = gprior.samplevar(ssq, M->nColUsed);
   }

   /* old code from ranf_cor for comparison
   void resid_correct(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] -= par[col] * matrixdata(coldata(obs),col);
   }
   
   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] += par[col] * matrixdata(coldata(obs),col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t rowlevel;
      for (size_t obs=0; obs < coldata.size(); obs++) {
         rowlevel = coldata(obs);
         rhs += matrixdata(rowlevel,col) * residPrec[obs] * resid[obs];
         lhs += matrixdata(rowlevel,col) * matrixdata(rowlevel,col) * residPrec[obs];
      }
   }

   */
   dataMatrix *M;
   dataFactor *F;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   
};

#endif /* modelMatrix_h */
