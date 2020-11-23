//
//  BayzR -- modelResp.h
//  model class for modelling response and residual variance. The base class is for the
//  standard continuous data case, there are derived classes for non-linear models.
// - the coldata is NumericVector (double)
// - hpar vector is size 1 to hold residual variance
// note: modelBase constructor has set parName to name of response variable
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelResp_h
#define modelResp_h

#include <Rcpp.h>
#include "modelBase.h"

class modelResp : public modelBase {
   
public:
   
   modelResp(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
            : modelBase(modelTerm, d, e, resp)
   {
      if(varColIndex[0] >= 0)
         Ydata = d[varColIndex[0]];
      else {  // this can be modified to use vector from the environment
         throw generalRbayzError("Not yet ready to use response ("+parName+") from R environment");
      }
      hpar.resize(1,1.0);
      hparName = "var.resid";
      // initialize resid and residPrec vectors; resid starts as the observed responses
      for (size_t obs=0; obs<Nresid; obs++) {
         resid[obs] = Ydata[obs];
         residPrec[obs] = 1.0;
      }
   }
   
   ~modelResp() {
   }

   void sample() {
      // Continuous data: sample() is only updating residual variance
      size_t obs;
      double sum=0.0;
      for (obs=0; obs<Nresid; obs++)
         sum += resid[obs]*resid[obs]*residPrec[obs];
      sum *= hpar[0];  // the sum was computed divided by the old variance!
      hpar[0] = gprior.samplevar(sum, Nresid);
      // the residPrec vector should be re-filled
      double temp = 1.0/hpar[0];
      for (obs=0; obs<Nresid; obs++)
         residPrec[obs] = temp;
   }
   
protected:
   Rcpp::NumericVector Ydata;

};

#endif /* modelResp_h */
