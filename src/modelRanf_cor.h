//
//  BayzR --- modelRanf_cor.hpp
//
// Computational class to to model ranf with correlation matrix (with V=K setting).
// This is a child-class of modelMatrix (change to Kernel?) because it works on a real matrix with regressions.
// - eigenvector data is already stored as coldata through constructor of modelTerm_realmat
// - eigenvalue data still needs to be set
// - hpar is one variance
// - par is set in modelTerm_realmat constructor and is size of coldata columns

//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_cor_h
#define modelRanf_cor_h

#include <Rcpp.h>
#include "modelMatrix.h"  // maybe should become modelKernel
#include "dataKernel.h"

class modelRanf_cor : public modelMatrix {  // maybe should become modelKernel

public:

   modelRanf_cor(Rcpp::DataFrame &d, size_t col) : modelMatrix(d,col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3];
      hparName = "var." + parName;
   }

   ~modelRanf_cor() {
   }
   

   void sample() {
      for(size_t k=0; k<nEvalUsed; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         lhs = lhs + (1.0/(eval[k]*hpar[0]));  // lhs with variance added
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<nEvalUsed; k++)
         ssq += par[k]*par[k]/eval[k];
      hpar[0] = gprior.samplevar(ssq,nEvalUsed);
   }
   
protected:

private:

   size_t nEvalUsed;  // this is in dataKernel

   Rcpp::NumericVector eval;
   Rcpp::IntegerVector update;
   
};

#endif /* modelRanf_cor_h */
