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

// update: ranfcor now also derive from modelFactor?
// the matrix data here was eigenvector data, it needs to come from a variance model now.

class modelRanf_cor : public modelFactor {

public:

   modelRanf_cor(dcModelTerm & modeldescr, modelBase * rmod)
         : modelFactor(modeldescr, rmod), regcoeff()
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
      if (modeldescr.varianceNames.size()==2) {  // combine with a second kernel if present
         kernelMatrix* K2 = new kernelMatrix(modeldescr.varianceObjects[1], modeldescr.varianceNames[1]);
         K.addKernel(K2);
         delete K2;
      }
      if (modeldescr.varianceNames.size()>2) {  // need to think if I can keep combining kernels with addKernel()
         throw(generalRbayzError("Not yet ready to combine more than 2 kernels for interaction"));
      }
      // note: modelled parameters are regressions on the eigenvectors and are stored in separate
      // vector regcoeff (size K->ncol). The par-vector has backtransfor to random effects
      // (eigen-vectors x regressions, size K->nrow).
      par.initWith(K->nrow,0.0l);
      parLabels = K->rownames; // I believe this makes a deep copy, but a shallow copy would be enough
      regcoeff.initWith(K->ncol, 0.0l);
      builObsIndex(obsIndex,F,K);
   }

   ~modelRanf_cor() {
   }
   
   void sample() {
/*    has to be revised completely, needs to get matrix part from a variance model
      update_regressions(TRUE, hpar[0]);
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]/weights[k];
      hpar[0] = gprior.samplevar(ssq, M->ncol);
*/
   }

   // prepForOutput puts the transform to breeding values in the par-vector
   void prepForOutput() {
      for(size_t row=0; row< M->nrow; row++) {
         par[row]=0.0l;
         for(size_t col=0; col<M->ncol; col++) {
            par[row] += M->data[col][row] * par[col];
         }
      }
   };

//   correlVarStr* varmodel;
   kernelMatrix* K;
   simpleDblVector regcoeff;
   std::vector<size_t> obsIndex;

};

#endif /* modelRanf_cor_h */
