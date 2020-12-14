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
   
   modelResp(std::string modelTerm, Rcpp::DataFrame &d, modelBase * rmod)
            : modelBase(modelTerm, d, rmod), weights()
   {
      // For now there is no check on the response vector to be correctly numerical, in future
      // rbayz main may need to make a triage for different types of response and then construct
      // appropriate response objects.
      Ydata = Rcpp::as<Rcpp::NumericVector>(varObjects[0]);
      Nresid = Ydata.size();
      par.initWith(Ydata);
      weights.initWith(Nresid,1.0l);
      resid = par.data;          // resid and residPrec are now alias for the par and weights
      residPred = weights.data;  // vectors inside a response object
      hpar.initWith(1,1.0l);
      // no parNames (labels for residuals) yet!
      hparName = "var.resid";
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

   simpleDblVector weights;
   double *resid, *residPrec;
   size_t Nresid;

protected:
   Rcpp::NumericVector Ydata;

};

#endif /* modelResp_h */
