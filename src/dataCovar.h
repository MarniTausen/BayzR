//
//  BayzR -- dataCovar.h
//  Storage class to hold a covariate
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef dataCovar_h
#define dataCovar_h

#include <Rcpp.h>

class dataCovar {
   
public:
   
   dataCovar(Rcpp::DataFrame &d, size_t col) {
      data = Rcpp::as<Rcpp::NumericVector>(d[col]);
   }
   
   ~dataCovar() {
   }
   
   Rcpp::NumericVector data;
   
};

#endif /* dataCovar_h */
