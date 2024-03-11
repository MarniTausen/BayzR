//  R/bayz
//  dataFactor.h - storing one or multiple interacting factors.
//  Interactions are recoded into a single factor with labels A1:B1:C1 etc., this is the setup
//  that can be used in fixed effects and uncorrelated random effects.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef dataFactor_h
#define dataFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "simpleVector.h"

class dataFactor {
public:
   dataFactor();
   dataFactor(Rcpp::DataFrame &d, size_t col);
   dataFactor(std::string varname);
   dataFactor(Rcpp::RObject Rcol);
   void setupFirstVariable(Rcpp::RObject col);
   void addVariable(Rcpp::DataFrame &d, size_t col);
   void addVariable(std::string varname);
   void addVariable(Rcpp::RObject Rcol);
   ~dataFactor() { }
   simpleIntVector data;
   std::vector<std::string> labels;
   int Nvar;  // The number of variables (interactions) in this factor
};

#endif /* dataFactor_h */
