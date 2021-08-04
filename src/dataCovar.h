//
//  BayzR -- dataCovar.h
//  Storage class to hold a covariate
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef dataCovar_h
#define dataCovar_h

#include <Rcpp.h>
#include "simpleVector.h"

class dataCovar : public simpleDblVector {

public:
   
   dataCovar(Rcpp::RObject x) : simpleDblVector() {
      Rcpp::NumericVector temp = Rcpp::as<Rcpp::NumericVector>(x);
      Rcpp::LogicalVector missing = Rcpp::is_na(temp);
      initWith(temp);
      size_t countnoNA=0;
      double sumnoNA=0.0l;
      for(size_t row=0; row<temp.size(); row++) {
         if(!missing[row]) {
            sumnoNA += data[row];
            countnoNA++;
         }
      }
      offset = sumnoNA / double(countnoNA);
      for(size_t row=0; row<temp.size(); row++) {
         if(missing[row])
            data[row] = 0.0;
         else
            data[row] -= offset;
      }
   }
   
   ~dataCovar() {
   }

   // when fitting a regression on covar, there will be a part offset * regcoeff missing from
   // the model intercept because that part has been substracted in the centering of the covar.
   double offset;
   
};

#endif /* dataCovar_h */
