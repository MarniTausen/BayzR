//
//  BayzR -- modelFixf.h
//  model class for fixed factor fit, derives from modelFactor.
//   - implements sample() method (see modelClasses.cpp for code)
//   - hpar is not used in fixf model
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelFixf_h
#define modelFixf_h

#include <Rcpp.h>
#include <cmath>
#include "modelFactor.h"

class modelFixf : public modelFactor {
   
public:
   
   modelFixf(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod) {
   }

   ~modelFixf() {
   }

   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=1; k<par.nelem; k++) {  // in fixf par[0] remains zero!, this runs from k=1
         if (lhs[k]>0)                     // if lhs is zero estimate will be set to 0
            par->val[k] = R::rnorm( (rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
         else
            par[k]=0.0;
      }
      resid_correct();
   }

};

#endif /* modelFixf_h */
