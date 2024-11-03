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
#include "parsedModelTerm.h"

#ifndef parVector_h
#define parVector_h

class parVector {

public:

   simpleDblVector Values;
   std::string Name="";
   std::string variables="";
   std::vector<std::string> Labels;
   std::string modelFunction="";
   std::string varianceStruct="";
   int traced=0;
   size_t nelem=0;
   double* val=0;  // convenience shortcut to retrieve parameter values as par->val[k]
   simpleDblVector postMean;
   simpleDblVector postVar;
   size_t count_collect_stats=0;
   parVector(parsedModelTerm & modeldescr, double initval);
   parVector(parsedModelTerm & modeldescr, double initval, std::string namePrefix);
   parVector(parsedModelTerm & modeldescr, double initval, Rcpp::CharacterVector& labels, std::string namePrefix);
   parVector(parsedModelTerm & modeldescr, double initval, Rcpp::CharacterVector& labels);
   parVector(parsedModelTerm & modeldescr, double initval, std::vector<std::string>& labels);
   void common_constructor_items(parsedModelTerm & modeldescr, std::string namePrefix);
   void collectStats();
   ~parVector() {   }
   
};

std::ostream& operator<<(std::ostream& os, const parVector& p);

#endif /* parVector_h */
