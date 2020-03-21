//
//  BayzR -- modelClasses.cpp
//
//  Created by Luc Janss on 24/01/2020.
//

#include <Rcpp.h>
#include <cmath>
#include "modelFactor.h"
#include "modelFixf.h"
#include "modelRanf.h"
#include "modelFreg.h"

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

void modelRanf::sample() {
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

void modelFreg::resid_correct() {
   for (size_t obs=0; obs < C->data.size(); obs++)
      resid[obs] -= par[0] * C->data[obs];
}

void modelFreg::resid_decorrect() {
   for (size_t obs=0; obs < C->data.size(); obs++)
      resid[obs] += par[0] * C->data[obs];
}

void modelFreg::collect_lhs_rhs() {
   lhs = 0.0; rhs=0.0;
   for (size_t obs=0; obs < C->data.size(); obs++) {
      rhs += residPrec[obs] * resid[obs] * C->data[obs];
      lhs += C->data[obs] * residPrec[obs] * C->data[obs];
   }
}

void modelFreg::sample() {
   resid_decorrect();
   collect_lhs_rhs();
   par[0] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
   resid_correct();
}
