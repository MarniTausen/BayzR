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
#include "indepVarStr.h"

class modelResp : public modelBase {
   
public:
   
   modelResp(dcModelTerm & modeldescr, modelBase * rmod)
            : modelBase(modeldescr, rmod)
   {
      // For now there is no check on the response vector to be correctly numerical, in future
      // rbayz main may need to make a triage for different types of response and then construct
      // appropriate response objects.
      Ydata = Rcpp::as<Rcpp::NumericVector>(varObjects[0]);
      Nresid = Ydata.size();
      par.initWith(Ydata);
      fname = "rp";
      parName = "resid";
      // no labels for residuals yet!
   }
   
   ~modelResp() {
   }

   void sample() {
/*    all this should now be done in indepVar objects
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
*/
   }

   indepVarStr* varModel;

protected:
   Rcpp::NumericVector Ydata;

};

#endif /* modelResp_h */
