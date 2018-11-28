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
#include "modelTerm_factor.h"
#include "priorClasses.h"

// random-correlated model term, built from ranf() model term with cor=K setting.
//

class modelTerm_ran_cor : public modelTerm_factor {

public:

   modelTerm_ran_cor(Rcpp::DataFrame &d, size_t col)  : modelTerm_factor(d, col) {
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
