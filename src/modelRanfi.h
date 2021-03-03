//
//  BayzR -- modelRanfi.h
//  model class for RAndom Factor fit with Indep variance structures, derives from modelFactor.
//   - implements sample() method (see modelClasses.cpp for code)
//   - hpar is set as "var..." of size 1
//   - sample() also includes updating variance
//  Note: random effect with correation is quite different code, it derives from modelMatrix.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanfi_h
#define modelRanfi_h

#include <Rcpp.h>
#include "modelFactor.h"

class modelRanfi : public modelFactor {

public:

   modelRanfi(dcModelTerm & modeldescr, modelBase * rmod)
         : modelFactor(modeldescr, rmod) {
      hpar.initWith(1,1.0l);
      hparName = "var." + parName;
   }

   ~modelRanfi() {
   }

   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=0; k<par.nelem; k++) {
         // random effect: add 1/hpar[0] in lhs
         par[k] = R::rnorm( (rhs[k]/(lhs[k]+(1/hpar[0]))), sqrt(1.0/(lhs[k]+(1/hpar[0]))));
      }
      resid_correct();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<par.nelem; k++)
         ssq += par[k]*par[k];
      hpar[0] = gprior.samplevar(ssq,par.nelem);
   }

};

#endif /* modelRanfi_h */
