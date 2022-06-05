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
   std::string funcName;
   std::string allVariableNames;
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   int varianceType;
   std::string varianceLinMod;
   std::vector<std::string> varianceParams;
   std::vector<std::string> varianceNames;
   std::vector<Rcpp::RObject> varianceObjects;
   int hierarchType; // 0=no, 1=simplified form index/matrix, 2=genuine
   std::string priorModel;
   std::string hierarchModel;
   std::string logging;
};

#endif /* dcModelTerm */
