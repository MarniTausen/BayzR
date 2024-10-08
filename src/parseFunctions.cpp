//
//  parseFunctions
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "nameTools.h"

using Rsize_t = long long int;
extern std::vector<std::string> Messages;
extern bool needStop;

void removeSpaces(std::string &s) {
   size_t pos;
   while((pos=s.find(' '))!=std::string::npos) s.erase(pos,1);
}

// Find matching closing bracket in a string, starting at fromPos which should be on an
// opening bracket. If during search a new opening bracket it found, this recurses on itself
// to first find the matching closing bracket of that one, etc...
// Works both on square [ and round ( brackets, searching for the closing one of the
// same kind.
// Returns the position of the matching closing bracket, or string::npos if failed to find it.
size_t findClosingBrack(std::string &s, size_t fromPos) {
   char closingBrack;
   if (s[fromPos] == '[') closingBrack = ']';
   else if (s[fromPos] == '(' ) closingBrack = ')';
   else
      throw(generalRbayzError("Wrong use of findClosingBrack"));
   while(s[fromPos] != closingBrack && fromPos != std::string::npos) {
      if(s[fromPos] == '[' || s[fromPos] == '(')
         fromPos = findClosingBrack(s,fromPos);
      if (fromPos != std::string::npos ) fromPos++;
   }
   return fromPos;
}

std::vector<std::string> splitString(std::string text, std::string splitchar) {
   std::vector<std::string> parts;
   if (text=="") {
      return parts;
   }
   size_t split=0, nextsplit;
   while (split != std::string::npos) {
      // find next split position starting from current split position
      nextsplit = text.find_first_of(splitchar,split);
      if (nextsplit == std::string::npos) {
         parts.push_back(text.substr(split,std::string::npos));
         split = nextsplit;
      }
      else {
         parts.push_back(text.substr(split,(nextsplit-split)));
         split = nextsplit + 1;
      }
   }
   return parts;
}

// Transform R model formula to a string
// Deparse can make multiple strings! 
std::string convertFormula(Rcpp::Formula f) {
   Rcpp::Function deparse("deparse");
   Rcpp::CharacterVector formparts = deparse(f,Rcpp::Named("width.cutoff")=100);
   std::string s;
   for(Rsize_t i=0; i<formparts.size(); i++ ) {
      s += Rcpp::as<std::string>(formparts[i]);
   }
   return(s);
}

// Return R object type (type of variable):
// 1: Factor, 2: IntegerVector, 3: NumericVector, 4: CharacterVector, 5: LogicalVector
// 6: (Numeric/double or integer) Matrix, 7: DataFrame, 8: List, 9: all other
int getVariableType(Rcpp::RObject x) {
   if( (Rcpp::is<Rcpp::NumericVector>(x) || Rcpp::is<Rcpp::IntegerVector>(x)) && Rf_isMatrix(x))
      return(6);
   if(Rcpp::is<Rcpp::NumericVector>(x) && !Rf_isMatrix(x))
      return(3);
   if(Rcpp::is<Rcpp::IntegerVector>(x)){
      if(Rf_isFactor(x)) return(1);
      else return(2);
   }
   if(Rcpp::is<Rcpp::CharacterVector>(x)) return(4);
   if(Rcpp::is<Rcpp::LogicalVector>(x)) return(5);
   if(Rcpp::is<Rcpp::DataFrame>(x)) return(7);
   if(Rcpp::is<Rcpp::List>(x)) return(8);
   return(9);
}

// Search and retrieve a variable 'name' by searching in the dataframe 'd' or in the R
// environment, and return it as an RObject. If not found return R_NilValue.
Rcpp::RObject getVariableObject(Rcpp::DataFrame &d, std::string name) {
   bool foundInEnv;
   Rcpp::RObject tempObject;
   Rcpp::Environment Renv;
   int colnr = findDataColumn(d, name);  // returns -1 if not found
   if(colnr>=0) {                        // found in data frame
      tempObject = Rcpp::as<Rcpp::RObject>(d[colnr]);
   }
   else {                               // search in R environment
      foundInEnv = Renv.exists(name);
      if(foundInEnv) tempObject = Renv[name];
      else return R_NilValue;           // not found in data frame or environment!
   }
   return tempObject;
}

// [!] New replaced getModelLHSTerms and getModelRHSTerms
// Split model-terms: splits model formula in individual pieces (separated by ~ and +).
// The return is a vector with strings, 0: response, 1: a 'mn(0)' or 'mn(1)' for intercept,
// 2 and further: the model terms from RHS.
// The intercept is formed with "mn()" so that mn comes in 'funcName' of the paresed model-term,
// and that is used in the main function to create the right modelling object.
std::vector<std::string> splitModelTerms(std::string mf) {
   std::vector<std::string> modelTerms;
   std::vector<std::string> lhsrhs = splitString(mf,"~");
   if(lhsrhs.size() != 2)
      throw(generalRbayzError("Model-formula has no (or multiple?) '~'"));
   if(lhsrhs[0].size()==0)
      throw(generalRbayzError("Model-formula has no response term(s)"));
   modelTerms.push_back(lhsrhs[0]);
   if(lhsrhs[1].size()==0) {        // RHS part is empty, insert default intercept
      modelTerms.push_back("mn(1)");    // and ready to return
      return modelTerms;
   }
   std::vector<std::string> RHSparts = splitString(lhsrhs[1],"+");
   if (RHSparts[0]=="0" || RHSparts[0]=="1") {  // There is an intercept specified
      if( RHSparts[0]=="0" )
         modelTerms.push_back("mn(0)");
      else
         modelTerms.push_back("mn(1)");
      for(size_t i=1; i<RHSparts.size(); i++)   // now there is RHSparts left from i=1
         modelTerms.push_back(RHSparts[i]);
   }
   else {                                       // no intercept specified:
      modelTerms.push_back("mn(1)");            // add an intercept and all
      for(size_t i=0; i<RHSparts.size(); i++)   // RHSparts counting from i=0
         modelTerms.push_back(RHSparts[i]);
   }
   return modelTerms;
}

/*
std::vector<std::string> getModelLHSTerms(std::string mf) {
   // split on ~, it must create two pieces that are the response term and list of explanatory terms
   std::vector<std::string> lhsrhs = splitString(mf,"~");
   if(lhsrhs.size() != 2)
      throw(generalRbayzError("Model-formula has no ~ to separate response and explanatory terms"));
   if(lhsrhs[0].size()==0)
      throw(generalRbayzError("Model-formula has no response term(s)"));
   // split the first part on + to make list of the LHS terms
   std::vector<std::string> parts = splitString(lhsrhs[0],"+");
   return parts;
}

std::vector<std::string> getModelRHSTerms(std::string mf) {
   // split on ~, it must create two pieces that are the response term and list of explanatory terms
   std::vector<std::string> lhsrhs = splitString(mf,"~");
   // split the RHS string on + to make list of explanatory model terms
   std::vector<std::string> parts = splitString(lhsrhs[1],"+");
   if(lhsrhs[1].size()==0) {   // RHS string was empty, insert default intercept
      parts.insert(parts.begin(),"1");
      return parts;
   }
   // RHS string not empty: check if there is an explicit "0" or "1", if not add a "1" term
   if( ! (parts[0]=="0" || parts[0]=="1") )
      parts.insert(parts.begin(),"1");
   return parts;
}
*/

std::vector<std::string> parseColNames(Rcpp::DataFrame & d, size_t col) {
   std::vector<std::string> names(7,"");
   Rcpp::CharacterVector AllColNames = d.names();
   std::string colName = Rcpp::as<std::string>(AllColNames[col]);
   // remove all spaces from the column name to facilitate parsing
   size_t pos;
   while((pos=colName.find(' '))!=std::string::npos) colName.erase(pos,1);
   // if there is a closing parenthesis also remove it
   if ( (pos = colName.find(')')) != std::string::npos)
      colName.erase(pos,1);
   // Search for opening parenthesis of function. If there is no function in
   // this column name, it is only a variable name and we're ready.
   if ( (pos = colName.find('(')) == std::string::npos) {
      names[1]=colName;
      return names;
   }
   // When continuing here, we're dealing with a function and pos is standing
   // on the opening parenthesis. Items inside the function are comma-separated.
   names[0]=colName.substr(0,pos);   // name of the function
   // split all the rest on comma, starting from pos++ (after opening parenthesis)
   pos++;
   std::vector<std::string> arguments = splitString(colName.substr(pos,std::string::npos), ",");
   if(arguments.size()==0)        // no arguments inside the function, we are ready
      return names;
   // now we're dealing with at least one term inside the function, it must start with variable name(s)
   std::vector<std::string> varNames = splitString(arguments[0], ":");
   for(size_t i=0; i<varNames.size() && i<3; i++ )
      names[i+1] = varNames[i];
   if (varNames.size() > 3)
      throw(generalRbayzError("There is a model term with more then 3-way interaction"));
   if (names[0]=="ranf" || names[0]=="rf") {    // for ranf and rf search for V= argument
      for(size_t i=1; i<arguments.size(); i++)
         if(arguments[i].substr(0,2)=="V=") names[3]=arguments[i].substr(2,std::string::npos);
   }
   if (names[0]=="ran2f") {   // for ran2f search for V1= and V2= in arguments
      for(size_t i=1; i<arguments.size(); i++) {
         if(arguments[i].substr(0,3)=="V1=") names[3]=arguments[i].substr(3,std::string::npos);
         if(arguments[i].substr(0,3)=="V2=") names[4]=arguments[i].substr(3,std::string::npos);
      }
   }
   return names;
}

// I found using argument DataFrame d, column names of the data frame are changed, by
//  defining as reference &d it is not changed. Looks like data frame is not a reference,
//  but making a copy.
std::string getWrapName(std::string modelTerm) {
   size_t pos1;
   pos1 = modelTerm.find('(');
   return(modelTerm.substr(0,pos1));
}

// Retrieving different parts of a model-term, syntax is:
//  funcname(varnames,V=...,prior=....,....)
// Note: spaces are removed from the model-formula, varnames, V=, prior= etc end with
// comma or closing parenthesis.
std::string getFuncName(std::string modelTerm) {
   size_t pos1 = modelTerm.find('(');
   if (pos1==std::string::npos)
      return "";
   else
      return modelTerm.substr(0, pos1);
}

std::string getVarNames(std::string modelTerm) {
   size_t pos1,pos2;
   pos1 = modelTerm.find('(');
   if(pos1==std::string::npos) {   // if there is no opening parenthesis, use the
      return(modelTerm);           // whole modelTerm as name and ready to return
   }
   pos2 = modelTerm.find_first_of("),",pos1);
   std::string temp = modelTerm.substr(pos1+1,pos2-pos1-1);
   return(temp);
}

std::string getVarDescr(std::string modelTerm) {
   size_t pos1,pos2;
   pos1 = modelTerm.find("V=");
   if(pos1==std::string::npos) {
      return "";
   }
   pos2 = pos1+1;  // pos2 starts scanning from character after V= until it has
                   // landed on a closing ')' or ',' that finishes the V= term, while
                   // skipping any nested opening-closing [] or () that it may find
                   // on the way as well as all commas inside nested brackets.
   while( ! (modelTerm[pos2] == ')' || modelTerm[pos2] == ',' || pos2==std::string::npos) ) {
      if(modelTerm[pos2]=='[' || modelTerm[pos2]=='(')
         pos2 = findClosingBrack(modelTerm, pos2);
      if(pos2 != std::string::npos) pos2++;
   }
   if(pos2 == std::string::npos) {  // somehow failed to get a variance description
      std::string s = modelTerm.substr(0,pos1);
      throw generalRbayzError("Syntax error in variance specification at: "+s);
      return "";
   }
   std::string temp = modelTerm.substr(pos1+2,pos2-pos1-2);
   return temp;
}

std::string getPriorDescr(std::string modelTerm) {
   size_t pos1,pos2;
   pos1 = modelTerm.find("prior=");
   if(pos1==std::string::npos) {
      return "";
   }
   pos2 = modelTerm.find_first_of("),",pos1);
   std::string temp = modelTerm.substr(pos1+6,pos2-pos1-1);
   return temp;
}

// Generic version to find any option in the model term, for instance
// use text "prior=", and this will return what follows after prior=,
// or empty string if "prior=" is not in the modelTerm.
std::string getOptionText(std::string modelTerm, std::string text) {
   size_t pos1,pos2;
   pos1 = modelTerm.find(text);
   if(pos1==std::string::npos) {
      return "";
   }
   pos2 = modelTerm.find_first_of("),",pos1);
   std::string temp = modelTerm.substr(pos1+text.length(),pos2-pos1-1);
   return temp;
}

// Check if keywords in user-specified options (stored in a map of (key,value) pairs) matches allowed
// list of options. Errors are stored in Messages list, but higher level needs to throw exception.
// Returns 0 for success, 1 for failure.
int checkOptions(std::map<std::string,std::string> userOpts, std::string allowedOpts) {
   std::istringstream ss(allowedOpts);
   std::map<std::string, int> allowedOptsMap;
   std::string s;                    
   while(ss >> s) allowedOptsMap[s];            // make a map of allowed options to search in it, but maybe
   std::map<std::string, int>::iterator f;      // regular find in unsorted vector<string> is also OK here
   for(std::map<std::string,std::string>::iterator p = userOpts.begin(); p!= userOpts.end(); ++p) {
      if( (f=allowedOptsMap.find(p->first)) == allowedOptsMap.end()) {
         Messages.push_back(" - Option not recognized (misspelled or not used here): " + p->first);
         needStop=true;
      }
   }
   if(needStop) return 1; else return 0;
}
