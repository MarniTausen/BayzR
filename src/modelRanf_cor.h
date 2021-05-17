//  BayzR --- modelRanf_cor.hpp
//  Computational class to to model random effect with correlation matrix (from using rn(..., V=K)).
//  Now also works with interactions and multiple kernels, but code now accepting ony max two kernels. 
//  Derives from modelFactor which sets-up the (interacting) factor data, the computational methods
//  though are more like modelMatrix and therefore parts of code are now same between modelMatrix and
//  ranf_cor. It would need multiple and virtual inheritance to re-use the modelMatrix methods here too.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_cor_h
#define modelRanf_cor_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "kernelMatrix.h"
#include "priorClasses.h"

// update: ranfcor now also derive from modelFactor?
// the matrix data here was eigenvector data, it needs to come from a variance model now.

class modelRanf_cor : public modelFactor {

public:

   modelRanf_cor(dcModelTerm & modeldescr, modelBase * rmod)
         : modelFactor(modeldescr, rmod), regcoeff(), gprior(modeldescr.priorModel)
   {
      hpar.initWith(1,1.0l);
      hparName = "var." + parName;
      // For the moment all variance objects must be kernels
      for(size_t i=0; i<modeldescr.varianceObjects.size(); i++) {
         if (modeldescr.varianceObjects[i]==R_NilValue) {
            throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
         }
      }
      // Get the first kernel and then add (making kronecker products) with second etc., if available
      K = new kernelMatrix(modeldescr.varianceObjects[0], modeldescr.varianceNames[0]);
      Rcpp::Rcout << "1\n";
      if (modeldescr.varianceNames.size()==2) {  // combine with a second kernel if present
         Rcpp::Rcout << "2\n";
         kernelMatrix* K2 = new kernelMatrix(modeldescr.varianceObjects[1], modeldescr.varianceNames[1]);
         Rcpp::Rcout << "3\n";
         K->addKernel(K2);
         Rcpp::Rcout << "4\n";
         delete K2;
      }
      Rcpp::Rcout << "5\n";
      if (modeldescr.varianceNames.size()>2) {  // need to think if I can keep combining kernels with addKernel()
         throw(generalRbayzError("Not yet ready to combine more than 2 kernels for interaction"));
      }
      // note: par-vector and names are set-up in modelFactor and follow levels of the data factor.
      // Here add a vector regcoeff (size K->ncol) to hold the regresssion on eigenvectors.
      regcoeff.initWith(K->ncol, 0.0l);
      builObsIndex(obsIndex,F,K);
   }

   ~modelRanf_cor() {
   }
   
   void sample() {

      // Update regressions on the eigenvectors
      double lhs, rhs;
      size_t matrixrow;
      double* colptr;
      for(size_t col=0; col < K->ncol; col++) {
         colptr = K->data[col];
         // residual de-correction for this evec column
         for (size_t obs=0; obs < F->data.nelem; obs++)
            resid[obs] += regcoeff[col] * colptr[obsIndex[obs]];
         // Make the lhs and rhs and update this column regression
         lhs = 0.0l; rhs=0.0l;
         for (size_t obs=0; obs < F->data.nelem; obs++) {
            matrixrow = obsIndex[obs];
            rhs += colptr[matrixrow] * residPrec[obs] * resid[obs];
            lhs += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
         }
         lhs += (1.0l / K->weights[col]) ;
         regcoeff[col] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         // residual correction for this column with updated regression
         for (size_t obs=0; obs < F->data.nelem; obs++)
            resid[obs] -= regcoeff[col] * colptr[obsIndex[obs]];
      }

      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t col=0; col< K->ncol; col++)
         ssq += regcoeff[col]*regcoeff[col]/K->weights[col];
      hpar[0] = gprior.samplevar(ssq, K->ncol);
   }

   // prepForOutput puts the transform to random effects in the par-vector
   void prepForOutput() {
      for(size_t row=0; row< K->nrow; row++) {
         par[row]=0.0l;
         for(size_t col=0; col<K->ncol; col++) {
            par[row] += K->data[col][row] * par[col];
         }
      }
   };

//   correlVarStr* varmodel;
   kernelMatrix* K;
   simpleDblVector regcoeff;
   std::vector<size_t> obsIndex;
   GenericPrior gprior;

};

#endif /* modelRanf_cor_h */
