//
//  BayzR -- parVector.h
//  Storage of parameter-vectors (every modelling object has one) with name for the whole vector,
//  labels for the elements, various other info, functions to collect posterior statistics, ...
//
//  Created by Luc Janss on 13/09/2023.
//

#include <Rcpp.h>
#include <string>
#include <vector>
#include "simpleVector.h"

#ifndef parVector_h
#define parVector_h

class parVector {

public:

   simpleDblVector parVals;
   std::string parName="";
   std::vector<std::string> parLabels;
   std::string modelFunction="";
   std::string variableList="";
   std::string varianceStruct="";
   size_t nelem=0;
   double* val=0;  // convenience shortcut to retrieve parameter values as par->val[k]
   simpleDblVector postMean;
   simpleDblVector postVar;
   size_t count_collect_stats=0;
   parVector(std::string parname, double initval);
   parVector(std::string parname, Rcpp::NumericVector initval, Rcpp::CharacterVector& labels);
   parVector(std::string parname, double initval, Rcpp::CharacterVector& labels);
   parVector(std::string parname, double initval, std::vector<std::string>& labels);
   void collecStats();
   void setParInfo(std::string fname, std::string vars, std::string varstruct);
   
};

std::ostream& operator<<(std::ostream& os, const parVector& p);

#endif /* parVector_h */
