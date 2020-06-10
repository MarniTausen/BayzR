//
//  BayzR -- modelResp.h
//  model class for modelling response and residual variance. The base class is for the
//  standard continuous data case, there are derived classes for non-linear models.
// - the coldata is NumericVector (double)
// - hpar vector is size 1 to hold residual variance
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelResp_h
#define modelResp_h

#include <Rcpp.h>
#include "modelBase.h"

class modelResp : public modelBase {
   
public:
   
   modelResp(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      Ydata = d[col];
      hpar.resize(1,1.0);
      parName = "unknown";  // par is not used!
      hparName = "var.resid";
      for (size_t obs=0; obs<resid.size(); obs++)
         // initialize resid vector by copying data-vector in it
         resid[obs] = Ydata[obs];
   }
   
   ~modelResp() {
   }

   void sample() {
      // Continuous data: sample() is only updating residual variance
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
   
protected:
   Rcpp::NumericVector Ydata;

};

#endif /* modelResp_h */
