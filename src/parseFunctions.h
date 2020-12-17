//
//  parseFunctions
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#include <Rcpp.h>
#include <vector>
#include <string>

#ifndef parseFunctions_h
#define parseFunctions_h

void removeSpaces(std::string &s);
std::vector<std::string> splitString(std::string text, char splitchar);
std::string convertFormula(Rcpp::Formula f);
int getVariableType(Rcpp::RObject x);
Rcpp::RObject getVariableObject(Rcpp::DataFrame &d, std::string name);
std::vector<std::string> getModelLHSTerms(std::string mf);
std::vector<std::string> getModelRHSTerms(std::string mf);
std::vector<std::string> parseColNames(Rcpp::DataFrame & d, size_t col);
std::string getWrapName(std::string modelTerm);
std::string getFuncName(std::string modelTerm);
std::string getVarNames(std::string modelTerm);
std::string getVarDescr(std::string modelTerm);
std::string getPriorDescr(std::string modelTerm);

#endif /* parseFunctions_h */
