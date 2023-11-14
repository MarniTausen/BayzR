//
//  dataFactor.h
//  rbayz
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
   // 'empty' constructor, it can be filled with addVariable()
   dataFactor() : data(), Nvar{0} { }
   // Other constructors
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
