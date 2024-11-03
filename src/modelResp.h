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
#include "parVector.h"

class modelResp : public modelBase {
   
public:
   
   modelResp(parsedModelTerm & modeldescr)
            : modelBase(), Y()
   {
      if(modeldescr.variableNames.size()>1) {
         throw generalRbayzError("Multiple response variables (" + modeldescr.variableString +
                                    ") cannot be handled this way");
      }
      if( ! (modeldescr.variableTypes[0]==2 || modeldescr.variableTypes[0]==3 )) {
         throw generalRbayzError("Response variable (" + modeldescr.variableString + 
                                    ") is not an R integer or numerical vector");
      }
      // The vectors in the response model are initialized as:
      //  tempY(origdata)  Y.data par(fitval)   resid
      //         y           y        0         y  
      //        NA           0        0         0
      // Note: the parameter vector that comes in the output has the fitted values, this gives more
      // complete information than the residuals would do for NAs. Residuals are in a member vector.
      Rcpp::NumericVector tempY = Rcpp::as<Rcpp::NumericVector>(modeldescr.variableObjects[0]);
      Y.initWith(tempY);
      Rcpp::IntegerVector labels_int = Rcpp::seq_len(tempY.size());
      Rcpp::CharacterVector labels = Rcpp::as<Rcpp::CharacterVector>(labels_int);
      par = new parVector(modeldescr,0.0l,labels,"fitval");
      resid = new parVector(modeldescr, 0.0l, labels,"resid");
      missing = Rcpp::is_na(tempY);
      for(size_t row=0; row<par->nelem; row++) {
         if(missing[row]) {
            Y.data[row] = 0.0l;
            resid->val[row] = 0.0l;
         }
         else
            resid->val[row] = Y.data[row];
      }
      // [ToDo] Move out when implementing more variance structures in residuals
      varModel = new idenVarStr(modeldescr,resid);
   }
   
   ~modelResp() {
      delete par;
      delete resid;
      delete varModel;
   }

   // readjust residuals to current fitted values.
   void readjResid() {
      for(size_t row=0; row<par->nelem; row++) {
         if(missing[row]) resid->val[row] = 0.0l;
         else resid->val[row] = Y.data[row] - par->val[row];
      }
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
   }

   void sampleHpars() {
      varModel->sample();
   }

   indepVarStr* varModel;
   Rcpp::LogicalVector missing;
   parVector* resid;

   simpleDblVector Y;

};

#endif /* modelResp_h */
