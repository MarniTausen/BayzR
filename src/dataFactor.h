//  R/bayz
//  dataFactor.h - storing one or multiple interacting factors with coding of the interaction levels
//     and making merged labels like "A1.B1.C1". 
//     Makes use of simpleFactor to code and store the single factors, the single factors are also kept
//     because some models make use of that (but for others it could be dropped ...).
//     If there is only one factor, dataFactor is basically a wrapper around one simpleFactor.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef dataFactor_h
#define dataFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "simpleFactor.h"

class dataFactor {
public:
   // the following constructors are only useful if there is only one factor in the model-term:
   dataFactor(Rcpp::DataFrame &d, size_t col);
   dataFactor(std::string varname);
   dataFactor(Rcpp::RObject Rcol);
   // for multiple interacting factors in the model-term (A:B:C) use this one:
   dataFactor(std::vector<Rcpp::RObject> variableObjects, std::vector<std::string> variableNames);
//   void setupFirstVariable(Rcpp::RObject col);
//   void addVariable(Rcpp::DataFrame &d, size_t col);
//   void addVariable(std::string varname);
//   void addVariable(Rcpp::RObject Rcol);
   ~dataFactor();
   simpleIntVector levcode;  // coded level info
   std::vector<std::string> labels;
   int Nvar;  // The number of variables (interactions) in this factor
   std::vector<simpleFactor *> factorList;
};

#endif /* dataFactor_h */
