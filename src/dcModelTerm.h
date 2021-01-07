//
//  dcModelTerm.h
//  A class to describe a DeConstructed Model-Term - where all elements of a model-term
//  are taken apart and stored in separate slots.
//
//  Created by Luc Janss on 29/12/2020.

#ifndef dcModelTerm_h
#define dcModelTerm_h

#include <Rcpp.h>
#include <vector>
#include <string>

class dcModelTerm {
public:
   dcModelTerm(std::string modelTerm, Rcpp::DataFrame &d);
   ~dcModelTerm() { }
   std::string funcName, allVariableNames;
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   int varianceType; // 0=not given, 1=with matrices, 2=with linear model
   int hierarchType; // 0=no, 1=simplified form with index variable, 2=genuine
   std::vector<std::string> varianceNames;  // list of matrices (names) for varianceType 1
   std::string varianceModel;  // the variance linear model description (varianceType 2)
   std::string priorModel;
   std::string hierarchModel;
};

#endif /* dcModelTerm */
