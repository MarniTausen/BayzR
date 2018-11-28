//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_random_h
#define modelTerm_random_h

#include <Rcpp.h>
#include "modelTerm_factor.h"
#include "priorClasses.h"

// Model-term for 'simple' random effect, this derives from modelTerm_factor, which has already prepared
// a lot of common things for model-terms handling a factor. 'Simple' random effect is IID variance,
// correlated random effects use ran-cor model-term.
// To add for ranf:
// - set hpar vector at size 1 to hold variance (more complex random effects could have >1 variance....)

class modelTerm_random : public modelTerm_factor {

public:

   modelTerm_random(Rcpp::DataFrame &d, size_t col)  : modelTerm_factor(d, col) {
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelTerm_random() {
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

#endif /* modelTerm_random_h */
