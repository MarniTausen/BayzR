//
//  BayzR --- modelMean.h
//
// Computational class for fitting mean (intercept)
// - this one has no data!
// - par is size 1 to hold mean
// - hpar is not used
// Note: the mean-term goes through the parsing because there is a "1" inserted in the
// model-terms, but it means the name set in the base-class constructor is "1" and needs
// to be changed.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMean_h
#define modelMean_h

#include <Rcpp.h>
#include "modelBase.h"


class modelMean : public modelBase {
   
public:

   modelMean(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelBase(modelTerm, d, e, resp) {
      par.resize(1,0);
      parName = "mean";
   }

   ~modelMean() {
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

#endif /* modelMean_h */
