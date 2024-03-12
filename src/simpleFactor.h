//  R/Bayz
//  RcppFactor.h - there is no 'native' factor in Rcpp, so I created one.
//
//  Created by Luc Janss on 11/03/2024.
//

#ifndef RcppFactor_h
#define RcppFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "simpleVector.h"

class RcppFactor {
public:
   RcppFactor(Rcpp::RObject col);
   ~RcppFactor();
   SimpleIntVector level;
   std::vector<std::string> labels;
};

#endif /* RcppFactor_h */
