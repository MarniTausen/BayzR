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
#include "modelMatrix.h"
#include "dataMatrix.h"

class modelRanf_cor : public modelMatrix {

public:

   modelRanf_cor(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelMatrix(modelTerm, d, e, resp)
   {
      hpar.resize(1,1);
      // modified: the parName is still the name(s) of variables, it needs some
      // modifications if multiple objects have the same variables. In old code
      // matrix names were added, but I think it gets too much to keep that system.
      hparName = "var." + parName;
      // note: modelMatrix has set-up the parameter vector to hold regressions on
      // the eigenvectors, but the output of this model class should be the random
      // effects on the original scale (eigen-vectors x regressions). So the par-vector
      // has size the number of matrix rows, and the prepForOutput() computes the
      // random effects from the regression coefficients when output is needed.
      par.resize(M->nrow,0);
      parLabels = M->labels; // I believe this makes a deep copy, but a shallow copy would be enough
   }

   ~modelRanf_cor() {
   }
   
   void sample() {
      update_regressions();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->ncol; k++)
         ssq += regcoeff[k]*regcoeff[k]/M->weights[k];
      hpar[0] = gprior.samplevar(ssq, M->ncol);
   }

   // prepForOutput puts the transform to breeding values in the par-vector
   void prepForOutput() {
      for(size_t row=0; row< M->nrow; row++) {
         par[row]=0.0l;
         for(size_t col=0; col<M->ncol; col++) {
            par[row] += M->data[col][row] * regcoeff[col];
         }
      }
   };
   
};

#endif /* modelRanf_cor_h */
