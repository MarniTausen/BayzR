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

   modelRanf_cor(Rcpp::DataFrame &d, size_t col) : modelMatrix(d,col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[4];
      hparName = "var." + parName;
      // note: modelMatrix has set-up the parameter vector to hold regressions on
      // the eigenvectors, but the output of this model class should be the random
      // effects on the original scale (eigen-vectors x regressions). So the par-vector
      // has size the number of matrix rows, and the prepForOutput() computes the
      // random effects from the regression coefficients when output is needed.
      par.resize(M->data.nrow(),0);
      parLabels = M->labels; // I believe this makes a deep copy, but a shallow copy would be enough
   }

   ~modelRanf_cor() {
   }
   
   void sample() {
      update_regressions();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->nColUsed; k++)
         ssq += regcoeff[k]*regcoeff[k]/M->weights[k];
      hpar[0] = gprior.samplevar(ssq, M->nColUsed);
   }

   // prepForOutput puts the transform to breeding values in the par-vector
   void prepForOutput() {
      for(size_t row=0; row< M->data.nrow(); row++) {
         par[row]=0.0l;
         for(size_t col=0; col<M->nColUsed; col++) {
            par[row] += M->data(row,col) * regcoeff[col];
         }
      }
   };
   
};

#endif /* modelRanf_cor_h */
