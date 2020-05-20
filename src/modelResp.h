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
         resid[obs] = Ydata[obs];  // initialize resid vector by copying data-vector in it
   }
   
   ~modelResp() {
   }

   void sample();
   
protected:
   Rcpp::NumericVector Ydata;

};

#endif /* modelResp_h */
