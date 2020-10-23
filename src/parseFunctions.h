//
//  parseColNames
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#include <Rcpp.h>
#include <vector>
#include <string>
#include "rbayzExceptions.h"

#ifndef parseColNames_h
#define parseColNames_h

// Function parseColNames retrieves names from a model-term: the function name, up to 3 variable names
// and up to 3 variance-matrices.
// Function getWrapName is a shortcut to return the function name
// Function getVarName is a shortcut to return the first variable
// Input is the model-frame, and the column number for which to extract the name (count from 0).
// Return is a vector of strings with 7 slots (0-6) for:
// 0: function name, 1,2,3: variable names, 4,5,6: variance matrix names.
// If anything is not available it is set "".

std::vector<std::string> splitString(std::string text, char splitchar) {
   std::vector<std::string> parts;
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


std::vector<std::string> parseModel(Rcpp::String mf) {
   std::string modelformula = Rcpp::as<std::string>(mf);
   // remove all spaces in the model formula
   size_t pos;
   while((pos=modelformula.find(' '))!=std::string::npos) modelformula.erase(pos,1);
   // split on ~, it must create two pieces that are the response term and list of explanatory terms
   std::vector<std::string> lhsrhs = splitString(modelformula,"~");
   if(lhsrhs.size() != 2)
      throw(generalRbayzError("Model-formula has no ~ to separate response and explanatory terms"));
   if(lhsrhs[0].size()==0)
      throw(generalRbayzError("Model-formula has no response term"));
   if(lhsrhs[1].size()==0)
      throw(generalRbayzError("Model-formula has no explanatory terms"));
   // split the RHS string on + to make list of explanatory model terms
   std::vector<std::string> parts = splitString(lhsrhs[1],"+");
   // Check for the first RHS term to be 0 or 1, and if not add a term '1' to set an intercept
   if(! (parts[0]=="0" || parts[0]=="1") )
      parts.push_front("1");
   // Finally put the response term also in front of the list of 'parts'
   parts.push_front(lhsrhs[0]);
   return parts;
}


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
   std::vector<std::string> arguments = splitString(colName.substr(pos,std::string::npos), ',');
   if(arguments.size()==0)        // no arguments inside the function, we are ready
      return names;
   // now we're dealing with at least one term inside the function, it must start with variable name(s)
   std::vector<std::string> varNames = splitString(arguments[0], ':');
   for(size_t i=0; i<varNames.size && i<3; i++ )
      names[i+1] = varNames[i];
   if (varNames.size>3)
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
std::string getWrapName(Rcpp::DataFrame & d, size_t col) {
   std::vector<std::string> names = parseColNames(d,col);
   return(names[0]);
}

std::string getVarName(Rcpp::DataFrame & d, size_t col) {
   std::vector<std::string> names = parseColNames(d,col);
   return(names[1]);
}


#endif /* parseColNames_h */
