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
#include <map>

class parsedModelTerm {
public:
   parsedModelTerm(std::string modelTerm, Rcpp::DataFrame &d);
   parsedModelTerm(std::string variableName, std::string VEdescr, Rcpp::DataFrame &d);
   ~parsedModelTerm() { }
   void parseModelTerm_step2(std::string funcName, std::string variables, 
                 std::string options, Rcpp::DataFrame &d);
   std::string funcName="";
   std::string shortModelTerm="";
   std::string variableString="";
   std::string variablePattern="";
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   std:map<std::string, std::string> user_options;  // could possibly replace varianceDescr and priormodDescr
   std::string varianceDescr="";
   std::string varianceStruct="";
   std::string varianceLinMod="";
   std::vector<std::string> varianceParams;
   std::vector<std::string> varianceNames;
   std::vector<Rcpp::RObject> varianceObjects;
   std::vector<std::string> varianceType;
   int hierarchType; // 0=no, 1=simplified form index/matrix, 2=genuine
   std::string priormodDescr="";
   std::string hierarchModel="";
   std::string logging="";
};

std::ostream& operator<<(std::ostream& os, const parsedModelTerm& p);

#endif /* parsedModelTerm */
