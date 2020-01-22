//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_mean_h
#define modelTerm_mean_h

#include <Rcpp.h>

// Model-term for fitting mean (intercept)
// - this one has no coldata!
// - par is size 1 to hold mean
// - hpar is not used

class modelTerm_mean : public modelTerm {
   
public:
   
   modelTerm_mean(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      par.resize(1,0);
      parName = "mean";
   }

   ~modelTerm_mean() {
   }

   void sample() {
      size_t obs, nobs=resid.size();
      double sum=0.0, temp=0.0;
      for (obs=0; obs<nobs; obs++) resid[obs] += par[0];
      for (obs=0; obs<nobs; obs++) {
         sum += resid[obs]*residPrec[obs];
         temp += residPrec[obs];
      }
      par[0] = R::rnorm((sum/temp), sqrt(1.0/temp));
      for (obs=0; obs<nobs; obs++) resid[obs] -= par[0];
   }
   
private:
};

#endif /* modelTerm_mean_h */
