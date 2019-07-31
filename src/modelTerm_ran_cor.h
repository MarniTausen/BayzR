//
//  modelTerm_ran_cor.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_ran_cor_h
#define modelTerm_ran_cor_h

#include <Rcpp.h>
#include "modelTerm_realmat.h"
#include "priorClasses.h"

// random-correlated model term, built from ranf() model term with cor=K setting.
// This is a child-class of modelTerm_realmat because it works on a real matrix with regressions.
// - in this model hpar is one variance
// - par is size number of columns of d (minus the ones with small eigenvalue?) to hold regression
//   coefficients

class modelTerm_ran_cor : public modelTerm_realmat {

public:

   modelTerm_ran_cor(Rcpp::DataFrame &d, size_t col, Rcpp::NumericMatrix &m)
                         : modelTerm_realmat(d, col, m) {
      hpar.resize(1,1);
      hparName = "var." + parName;
      
      // add setting of NumericVector and NumericMatrix that match eigen-values and vectors
      
      
   }

   ~modelTerm_ran_cor() {
   }

   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=0; k<par.size(); k++) {
         // random effect: add 1/hpar[0] in lhs
         par[k] = R::rnorm( (rhs[k]/(lhs[k]+(1/hpar[0]))), sqrt(1.0/(lhs[k]+(1/hpar[0]))));
      }
      resid_correct();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<par.size(); k++)
         ssq += par[k]*par[k];
      hpar[0] = gprior.samplevar(ssq,par.size());
   }

private:

};

#endif /* modelTerm_ran_cor_h */
