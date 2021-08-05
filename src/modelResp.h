//
//  BayzR -- modelResp.h
//  model class for modelling response and residual variance. The base class is for the
//  standard continuous data case, there are derived classes for non-linear models.
// note: modelBase constructor has set parName to name of response variable
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelResp_h
#define modelResp_h

#include <Rcpp.h>
#include "modelBase.h"
#include "indepVarStr.h"
#include "simpleVector.h"

class modelResp : public modelBase {
   
public:
   
   modelResp(dcModelTerm & modeldescr, modelBase * rmod)
            : modelBase(modeldescr, rmod), fit(), Y()
   {
      // For now there is no check on the response vector to be correctly numerical, in future
      // rbayz main may need to make a triage for different types of response and then construct
      // appropriate response objects.
      Rcpp::NumericVector tempY = Rcpp::as<Rcpp::NumericVector>(varObjects[0]);
      Y.initWith(tempY);
      par.initWith(tempY);
      missing = Rcpp::is_na(tempY);
      for(size_t row=0; row<par.nelem; row++) {
         if(missing[row]) {
            par.data[row] = 0.0l;
            Y.data[row] = 0.0l;
         }
      }
      fname = "rp";
      parName = "resid";
      fit.initWith(par.nelem, 0.0l);
      // no labels for residuals yet!
   }
   
   ~modelResp() {
   }

   void sample() {
      // Store fitted values and re-sample the Ydata / residuals for missing data
      for(size_t row=0; row<par.nelem; row++) {
         if(missing[row]) {
            fit.data[row] = Y.data[row] - par.data[row];
            par.data[row] = R::rnorm( 0.0l, sqrt(1.0/varModel->weights[row]));
            Y.data[row] = fit.data[row] + par.data[row];
         }
         else {
            fit.data[row] = Y.data[row] - par.data[row];
         }
      }

/*    variance update is now moded to indepVar objects
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
   Rcpp::LogicalVector missing;
   simpleDblVector fit;


protected:
   simpleDblVector Y;

};

#endif /* modelResp_h */
