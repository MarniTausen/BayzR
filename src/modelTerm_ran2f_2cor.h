//
//  modelTerm_ran2f_2cor.hpp
//  rbayz
//
//  Created by Luc Janss on 18/01/2020.
//

#ifndef modelTerm_ran2f_2cor_h
#define modelTerm_ran2f_2cor_h

#include <Rcpp.h>
#include "modelTerm_2factor.h"
#include "priorClasses.h"

// Model-term for interaction between two random factors both with covariance-structures.

class modelTerm_ran2f_2cor : public modelTerm_2factor {

public:

   modelTerm_ran2f_2cor(Rcpp::DataFrame &d, size_t col)  : modelTerm_2factor(d, col, col2) {
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelTerm_ran2f() {
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

   // this is code from ran2f that collects all lhs and rhs at once, probably
   // needs to be changed to collect one level at-a-time like ran_cor.
   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.size(); k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs<intdata.size(); obs++) {
         k=intdata[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

};

#endif /* modelTerm_ran2f_2cor_h */
