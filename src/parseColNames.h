//
//  parseColNames
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#include <Rcpp.h>
#include <vector>
#include <string>

#ifndef parseColNames_h
#define parseColNames_h

// Two functions to separate parts of the column names in the model-frame.
// The column names have syntax:
//   variable
//   wrapfunc(variable)
//   wrapfunc(variable, option)
//   wrapfunc( variable , option)   <- extra spaces may occurr around variable
// Function getWrapName extracts the 'wrapfunc' part (returns "" if there is no)
// Function getVarName extracts the 'variable' part.
// Input is the model-frame, and the column number for which to extract the name (count from 0).

//! I found using argument DataFrame d, column names of the data frame are changed, by
//  defining as reference &d it is not changed. Looks like data frame is not a reference,
//  but making a copy.
std::string getWrapName(Rcpp::DataFrame & d, size_t col) {
   Rcpp::CharacterVector colNames = d.names();
   std::string colName = Rcpp::as<std::string>(colNames[col]);
   size_t leftparenth = colName.find('(');
   std::string wrapName = (leftparenth == std::string::npos) ? "" : colName.substr(0,leftparenth);
   return(wrapName);
}

//! can add some checks, for instance rightpos should not get <= leftpos, and maybe more exceptions ...
std::string getVarName(Rcpp::DataFrame & d, size_t col) {
   Rcpp::CharacterVector colNames = d.names();
   std::string colName = Rcpp::as<std::string>(colNames[col]);
   // First, set leftpos for the variable name at 0 (if there is no left parenthesis), or 1 after left parenthesis
   size_t leftpos = (colName.find('(') == std::string::npos) ? 0 : colName.find('(') + 1 ;
   // To accommodate for a space before the variable name, increase leftpos as long as it is a space
   while(colName[leftpos]==' ' && leftpos < colName.length()) leftpos++;
   // And then find the right-position as the first space, comma or right parenthesis after leftpos
   size_t rightpos = colName.find_first_of(" ,)",leftpos);
   std::string varName = colName.substr(leftpos,(rightpos-leftpos));
   return(varName);
}

#endif /* parseColNames_h */
