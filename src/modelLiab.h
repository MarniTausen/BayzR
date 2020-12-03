//
//  BayzR -- modelLiab.h
//  model class for modelling liabilities for categorical response data- a derived class of modelResp.
//  - in the parent class modelResp the input column is stored as continuous data in Ydata, but
//    here the same input is again stored as integer data in CatData. The Ydata will be used for
//    storing liabilies.
// - par vector is used for the thresholds excluding the -INF and +INF thresholds
// - hpar vector is size 1 to hold residual variance
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelLiab_h
#define modelLiab_h

#include <Rcpp.h>
#include <string>
#include <utility>
#include "modelResp.h"

class modelLiab : public modelResp {
   
public:
   
   modelLiab(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelResp(modelTerm, d, e, resp) {
      catData = Rcpp::as<Rcpp::IntegerVector>(d[varColIndex[0]]);
      Rcpp::IntegerVector categories = Rcpp::sort_unique(catData);
      par.initWith(categories.size()-1, 0.0l);
      parName = "thresh";
      std::string s;  // Make threshold names
      for (size_t cat=0; cat<par.size(); cat++) {
         s = std::to_string(categories[cat]) + ":" + std::to_string(categories[cat+1]);
         parLabels.push_back(s);
      }
      // there is a second set of thresholds including the -INF and +INF ones for internal computations
      thresholds.resize(categories.size()+1);
      for(size_t t=0; t<=thresholds.size(); t++) thresholds[t]=double(t-1);
      // Fill Ydata and residuals with a first sample of liabilities
      std::pair<double, double> bounds;
      for (size_t obs=0; obs<resid.size(); obs++) {
         bounds = liabBounds(catData[obs]);  // obtain liability boundaries for this obs, this
                                             // also depends on current liab or residual?
         Ydata[obs] = R::runif(bounds.first,bounds.second);
         resid[obs] = Ydata[obs];
      }
   }
   
   ~modelLiab() {
   }

   std::pair<double, double> liabBounds(int i) {
      std::pair<double, double> bounds;
      // not yet finished here ...
      return bounds;
   }

   void sample();

private:
   Rcpp::IntegerVector catData;
   std::vector<double> thresholds;
   std::pair<double, double> liabBounds(int i);

};

#endif /* modelLiab_h */
