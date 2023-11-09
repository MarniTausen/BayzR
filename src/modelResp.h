//
//  BayzR -- modelResp.h
//  model class for modelling response and residual variance. The base class is for the
//  standard continuous data case, there are derived classes for non-linear models.
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
   
   modelResp(parsedModelTerm & modeldescr)
            : modelBase(), fit(), Y()
   {
      if(resp.variableNames.size()>1) {
         throw generalRbayzError("Multiple response variables (" + resp.variableString + ") cannot be handled this way");
      }
      if( ! (resp.variableTypes[0]==2 || resp.variableTypes[0]==3 )) {
         throw generalRbayzError("Response variable (" + resp.variableString + ") is not an R integer or numerical vector");
      }
      Rcpp::NumericVector tempY = Rcpp::as<Rcpp::NumericVector>(resp.variableObjects[0]);
      Y.initWith(tempY);
      Rcpp::IntegerVector labels_int = Rcpp::seq_len(tempY.size());
      Rcpp::CharacterVector labels = Rcpp::as<Rcpp::CharacterVector>(labels_int);
      par=new parVector(modeldescr,tempY,labels);
//    overwrite default naming made in modelBase constructor
      par->parName="fval."+modeldescr.variableString;
      missing = Rcpp::is_na(tempY);
      for(size_t row=0; row<par->nelem; row++) {
         if(missing[row]) {
            par->val[row] = 0.0l;
            Y.data[row] = 0.0l;
         }
      }
      resid.initWith(par->nelem, 0.0l);
      varModel = new idenVarStr(resp,this->par);
   }
   
   ~modelResp() {
   }

   void sample() {
      // Store fitted values and re-sample the Ydata and residuals for missing data
      // Fit, resid and Y are all needed because other objects modify the residuals, and with fit and Y
      // it is possible to keep track of these modifications and collect fitted values for missing data.
      for(size_t row=0; row<par->nelem; row++) {
         if(missing[row]) {
            par->val[row] = Y.data[row] - resid->val[row];
            resid->val[row] = R::rnorm( 0.0l, sqrt(1.0/varModel->weights[row]));
            Y.data[row] = par->val[row] + resid->val[row];
         }
         else {
            par->val[row] = Y.data[row] - resid->val[row];
         }
      }
/*    variance update is now moved to indepVar objects
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
      varModel->sample();
   }

   indepVarStr* varModel;
   Rcpp::LogicalVector missing;
   simpleDblVector resid;


protected:
   simpleDblVector Y;

};

#endif /* modelResp_h */
