//
//  BayzR -- modelClasses.cpp
//
//  Created by Luc Janss on 24/01/2020.
//

#include <Rcpp.h>
#include <cmath>
#include "modelFactor.h"
#include "modelFixf.h"

void modelFactor::resid_correct() {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] -= par[F->data[obs]];
}

void modelFactor::resid_decorrect() {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] += par[F->data[obs]];
}

void modelFactor::collect_lhs_rhs() {
   size_t k;
   for(k=0; k<par.size(); k++) {
      rhs[k] = 0.0;
      lhs[k] = 0.0;
   }
   for (size_t obs=0; obs < F->data.size(); obs++) {
      k=F->data[obs];
      rhs[k] += residPrec[obs] * resid[obs];
      lhs[k] += residPrec[obs];
   }
}

void modelFixf::sample() {
   resid_decorrect();
   collect_lhs_rhs();
   for(size_t k=1; k<par.size(); k++) {  // in fixf par[0] remains zero!, this runs from k=1
      if (lhs[k]>0)                      // if lhs (inverse variance) is zero estimate will be set to 0
         par[k] = R::rnorm( (rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
      else
         par[k]=0.0;
   }
   resid_correct();
}

