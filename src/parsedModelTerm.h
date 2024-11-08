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
   std::map<std::string, std::string> options; 
   std::string varianceStruct="";
   std::string varianceLinMod="";
   // vectors from variance-term parsing
   std::vector<std::string> varName;
   std::vector<std::string> varOption;
   std::vector<std::string> varVariable;
   std::vector<Rcpp::RObject> varObject;
   std::vector<std::string> varType;
   int hierarchType; // 0=no, 1=simplified form index/matrix, 2=genuine
   std::string hierarchModel="";
};

std::ostream& operator<<(std::ostream& os, parsedModelTerm& p);

#endif /* parsedModelTerm */
