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

// random-correlated model term, built from ranf() model term with V=K setting.
// This is a child-class of modelTerm_realmat because it works on a real matrix with regressions.
// - eigenvector data is already stored as coldata through constructor of modelTerm_realmat
// - eigenvalue data still needs to be set
// - hpar is one variance
// - par is set in modelTerm_realmat constructor and is size of coldata columns

class modelTerm_ran_cor : public modelTerm_realmat {

public:

   modelTerm_ran_cor(Rcpp::DataFrame &d, size_t col) : modelTerm_realmat(d, col, 2) {
      hpar.resize(1,1);
      hparName = "var." + parName;
      Rcpp::RObject thiscol = d[col];
      eval = Rcpp::as<Rcpp::NumericVector>(thiscol.attr("evalues"));
      nPosEval=0;
      for (size_t i=0;  i<eval.size() && eval[i] >= 0.001; i++)
         nPosEval++;
      
   }

   ~modelTerm_ran_cor() {
   }

   void sample() {
      double temp;
      for(size_t k=0; k<nPosEval; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         temp = lhs + (1.0/(eval[k]*hpar[0]));  // lhs with variance added
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<nPosEval; k++)
         ssq += par[k]*par[k]/eval[k];
      hpar[0] = gprior.samplevar(ssq,nPosEval);
   }

private:

   Rcpp::NumericVector eval;
   Rcpp::IntegerVector update;
   size_t nPosEval;
   
};

#endif /* modelTerm_ran_cor_h */
