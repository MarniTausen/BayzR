//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_resp_h
#define modelTerm_resp_h

#include <Rcpp.h>
#include "parseColNames.h"

// Model-term for response column.
// The response column updates residuals and liabilities (not yet implemented) and residual variance.
// - the coldata is NumericVector (double)
// - par vector is now not used (size will remain zero), it could be used for liabilities later
// - hpar vector is size 1 to hold residual variance

class modelTerm_resp : public modelTerm {
   
public:
   
   modelTerm_resp(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      coldata = d[col];
      hpar.resize(1,1.0);
      parName = "resid";
      hparName = "var.resid";
      for (size_t obs=0; obs<resid.size(); obs++)
         resid[obs] = coldata[obs];  // initialize resid vector by copying data-vector in it
   }
   
   ~modelTerm_resp() {
   }

   void sample() {
      // Re-sample residuals and liabilities (not yet done ...) ... and residuals for missing data?
      // ....
      // Update residual variance
      size_t obs, nobs=resid.size();
      double sum=0.0;
      for (obs=0; obs<nobs; obs++)
         sum += resid[obs]*resid[obs]*residPrec[obs];
      sum *= hpar[0];  // the sum was computed divided by the old variance!
      hpar[0] = gprior.samplevar(sum, nobs);
      // the residPrec vector should be re-filled
      double temp = 1.0/hpar[0];
      for (obs=0; obs<nobs; obs++)
         residPrec[obs] = temp;
   }
   
private:
   Rcpp::NumericVector coldata;

};

#endif /* modelTerm_resp_h */
