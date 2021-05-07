//
//  BayzR --- modelRanf_cor.hpp
//
//  Computational class to to model ranf with correlation matrix (with V=K setting).
//  This used all methods from modelMatrix, only some parameter names need to be set.
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
         : modelFactor(modeldescr, rmod)
   {
      hpar.initWith(1,1.0l);
      hparName = "var." + parName;
      // For the moment all variance objects must be kernels
      for(size_t i=0; i<modeldescr.varianceObjects.size(); i++) {
         if (modeldescr.varianceObjects[i]==R_NilValue) {
            throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
         }
      }
      // Get the first kernel
      K = new kernelMatrix(modeldescr.varianceObjects[0], modeldescr.varianceNames[0]);
      if (modeldescr.varianceNames.size()==2) {  // combine with a second kernel if present
         kernelMatrix* K2 = new kernelMatrix(modeldescr.varianceObjects[1], modeldescr.varianceNames[1]);
         K.addKernel(K2);
         delete K2;
      }
      if (modeldescr.varianceNames.size()>2) {  // need to think if I can keep combining kernels with addKernel()
         throw(generalRbayzError("Not yet ready to combine more than 2 kernels for interaction"));
      }

      // note: modelMatrix has set-up the parameter vector to hold regressions on
      // the eigenvectors, but the output of this model class should be the random
      // effects on the original scale (eigen-vectors x regressions). So the par-vector
      // has size the number of matrix rows, and the prepForOutput() computes the
      // random effects from the regression coefficients when output is needed.
      par.initWith(M->nrow,0.0l);
      parLabels = M->rownames; // I believe this makes a deep copy, but a shallow copy would be enough
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

};

#endif /* modelRanf_cor_h */
