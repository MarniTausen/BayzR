//
//  BayzR -- modelRanf.h
//  model class for random factor fit without correlation, derives from modelFactor.
//   - implements sample() method (see modelClasses.cpp for code)
//   - hpar is set as "var..." of size 1
//   - sample() also includes updating variance
//  Note: random effect with correation is quite different code, it derives from modelMatrix.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_h
#define modelRanf_h

#include <Rcpp.h>
#include "modelFactor.h"

class modelRanf : public modelFactor {

public:

   modelRanf(Rcpp::DataFrame &d, size_t col)  : modelFactor(d, col) {
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelRanf() {
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

};

#endif /* modelRanf_h */
