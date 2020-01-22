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

// Function parseColNames retrieves names from a model-term: the function name, up to 2 variable names
// and matrix-names set at V= or V1= and V2=.
// Function getWrapName is a shortcut to return the function name
// Function getVarName is a shortcut to return the first variable
// Input is the model-frame, and the column number for which to extract the name (count from 0).
// Return is a vector of strings with 5 slots for:
// 0: function name, 1:1st variable name, 2:2nd variable name, 3: 1st V=name, 4: 2nd V=name.
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

std::vector<std::string> parseColNames(Rcpp::DataFrame & d, size_t col) {
   std::vector<std::string> names(5,"");
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
   else {  // there is a function, can store first part as wrapfunc name
      names[0]=colName.substr(0,pos);
   }
   // split all the rest on comma, starting from pos++ (after opening parenthesis)
   pos++;
   std::vector<std::string> arguments = splitString(colName.substr(pos,std::string::npos), ',');
   if (arguments.size() > 0)
      names[1]=arguments[0];       // the first function argument is always the (first) variable name
   if (arguments.size() > 1 && arguments[1].find('=') == std::string::npos ) {
      names[2]=arguments[1];   // if 2nd argument has no '=' it is probably a second variable name
   }
   if (names[0]=="ranf") {    // for ranf search for V= argument (it is not required to be there)
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
