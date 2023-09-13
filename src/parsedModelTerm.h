//  parsedModelTerm.h
//  A class to describe a parsed Model-Term - where all elements of a model-term
//  are taken apart and stored in separate slots.
//
//  Created by Luc Janss on 29/12/2020.

#ifndef parsedModelTerm_h
#define parsedModelTerm_h

#include <Rcpp.h>
#include <vector>
#include <string>

class parsedModelTerm {
public:
   parsedModelTerm(std::string modelTerm, Rcpp::DataFrame &d);
   ~parsedModelTerm() { }
   std::string funcName="";
   std::string shortModelTerm="";
   std::string variableString="";
   std::string variablePattern="";
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   std::string varianceDescr="";
   std::string varianceStruct;
   std::string varianceLinMod="";
   std::vector<std::string> varianceParams;
   std::vector<std::string> varianceNames;
   std::vector<Rcpp::RObject> varianceObjects;
   std::vector<int> varianceKernelType;
   int hierarchType; // 0=no, 1=simplified form index/matrix, 2=genuine
   std::string priormodDescr="";
   std::string hierarchModel="";
   std::string logging="";
};

#endif /* parsedModelTerm */
