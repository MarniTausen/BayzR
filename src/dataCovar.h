//
//  BayzR -- dataCovar.h
//  Storage class to hold a covariate
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef dataCovar_h
#define dataCovar_h

#include <Rcpp.h>
#include "simpleVector"

class dataCovar : public simpleDblVector {

public:
   
   dataCovar(Rcpp::RObject x) : simpleDblVector() {
      Rcpp::NumericVector temp = Rcpp::as<Rcpp::NumericVector>(x);
      initWith(temp);
   }
   
   ~dataCovar() {
   }
   
   Rcpp::NumericVector data;
   
};

#endif /* dataCovar_h */
