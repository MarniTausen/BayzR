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
#include "simpleVector.h"
#include "modelBase.h"
#include "nameTools.h"

class modelMatrix : public modelBase {
   
public:
   
   modelMatrix(dcModelTerm & modeldescr, modelBase * rmod)
         : modelBase(modeldescr, rmod)
   {
      // For now only allowing a matrix input where there is an index variable (model
      // made with id/matrix). It could be extended to allow for no id, so that matrix needs to
      // be aliged 1:1 with data records, then the 'id' is bascially a 1:1 link.
      if( ! (varType[0]==1 && varType[1]==6 && hierType==1) )
         throw generalRbayzError("rr() model not supported, can now only deal with <factor>/<matrix> input");
      F = new dataFactor(varObjects[0]);
      M = new dataMatrix(varObjects[1], varNames[1]);
      par.initWith(M->ncol,0.0l);
      weights.initWith(M->ncol,1.0l);
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
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] -= par[col] * colptr[obsIndex[obs]];
   }

   void resid_decorrect(size_t col) {
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] += par[col] * colptr[obsIndex[obs]];
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t matrixrow;
      double * colptr = M->data[col];
      for (size_t obs=0; obs < F->data.nelem; obs++) {
         matrixrow = obsIndex[obs];
         rhs += colptr[matrixrow] * residPrec[obs] * resid[obs];
         lhs += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
      }
   }

// Update regressions for three different variance models/specifications
   void update_regressions(bool useWeights, double var) {
      //  - only use weights (var is set <0)
      if (useWeights && var <= 0) {
         for(size_t k=0; k < M->ncol; k++) {
            resid_decorrect(k);
            collect_lhs_rhs(k);
            lhs = lhs + (1.0/( weights[k] ));
            par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));  // Note weights are not precisions here
            resid_correct(k);
         }
      }
      //  - using weights and an extra variance scaling as weights*var
      else if (useWeights && var >0) {
         for(size_t k=0; k < M->ncol; k++) {
            resid_decorrect(k);
            collect_lhs_rhs(k);
            lhs = lhs + (1.0/( weights[k]*var ));
            par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));  // Note weights are not precisions here
            resid_correct(k);
         }
      }
      //  - not using weights
      else if (!useWeights && var> 0) {
         for(size_t k=0; k < M->ncol; k++) {
            resid_decorrect(k);
            collect_lhs_rhs(k);
            lhs = lhs + (1.0/( var ));
            par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));  // Note weights are not precisions here
            resid_correct(k);
         }
      }
      else
         throw generalRbayzError("Invalid useWeights and var settings in modelMatrix::update_regressions");
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
   simpleDblVector weights;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   std::vector<size_t> obsIndex;

};

#endif /* modelMatrix_h */
