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
#include "nameTools.h"

class modelMatrix : public modelBase {
   
public:
   
   modelMatrix(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelBase(modelTerm, d, e, resp)
   {
      size_t col = varColIndex[0]; 
      Rcpp::RObject col_asRObject = d[col];
      M = new dataMatrix(col_asRObject);
      F = new dataFactor(d, col);
      regcoeff.resize(M->nColUsed,0);
      // here level-names for the regression coefficients are not filled in parLevelNames,
      // therefore automatically the levels 1,2,3,... are inserted in collectParInfo.
      // It would be nicer to collect the column names (if available), or else fill with "col"+1,2,3...
      builObsIndex(obsIndex,F,M);
      lhs = 0.0l;
      rhs = 0.0l;
   }
   
   ~modelMatrix() {
      delete M;
      delete F;
   }
   
   void resid_correct(size_t col) {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] -= regcoeff[col] * M->data(obsIndex[obs],col);
   }

   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] += regcoeff[col] * M->data(obsIndex[obs],col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t matrixrow;
      for (size_t obs=0; obs < F->data.size(); obs++) {
         matrixrow = obsIndex[obs];
         rhs += M->data(matrixrow,col) * residPrec[obs] * resid[obs];
         lhs += M->data(matrixrow,col) * M->data(matrixrow,col) * residPrec[obs];
      }
   }

   void update_regressions() {
      for(size_t k=0; k < M->nColUsed; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         lhs = lhs + (1.0/( M->weights[k]*hpar[0] ));   // lhs with variance added
         regcoeff[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));  // Note weights still stored as variances, not inverse
         resid_correct(k);
      }
   }
   
   // Here no sample() yet, modelMatrix remains virtual. The derived classes implement sample()
   // by combining update_regressions() with update of hyper-paramters for that derived class.

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
   std::vector<double> regcoeff;
   std::vector<size_t> obsIndex;

};

#endif /* modelMatrix_h */
