//  dcModelTerm.cpp

#include <Rcpp.h>
#include <vector>
#include <string>
#include "dcModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

dcModelTerm::dcModelTerm(std::string modelTerm, Rcpp::DataFrame &d) :
   varianceType(0), hierarchType(0), varianceModel(""), priorModel(""), hierarchModel("")
{
   // funcName will be empty if there is no opening parenthesis
   size_t pos;
   if ( (pos = modelTerm.find('(')) == std::string::npos)
      funcName = "";
   else
      funcName = modelTerm.substr(0, pos);
   // Variables name(s) need some pre-treatment to handle hierarchy specifications
   std::string tempnames = getVarNames(modelTerm);
   pos = tempnames.find('/');             // search '/' for hierarchical models
   if( pos != std::string::npos) {
      size_t pos2 = tempnames.find('(');  // proper hierarchy needs parentheses somewhere,
      if (pos2 == std::string::npos) {    // without any parenthesis it is evaluated as
         hierarchType=1;                  // the simplified form variable/variable and
         tempnames[pos]=':';              // the '/' is replaced by ':' for further processing
      }
      else {                              // This is proper hierarchy, the hierarchical part
         hierarchType=2;                  // is removed and placed in hierarchModel string
         hierarchModel =                  // for later processing
                 tempnames.substr(pos+1, std::string::npos);
         tempnames.erase(pos, std::string::npos);
      }
   }
   // Here the variable names get in their final storage as variable names in a vector
   allVariableNames = tempnames;
   variableNames = splitString(tempnames,':');
   // For every variable get an RObject pointing to it (whether it is from the data frame
   // or from R environment), and also store the type (factor, numeric, etc.) of the variable.
   for(size_t i=0; i<variableNames.size(); i++) {
      if (variableNames[i]=="1" || variableNames[i]=="0") {
         variableObjects.push_back(R_NilValue);
         variableTypes.push_back(0);
      }
      else {
         variableObjects.push_back(getVariableObject(d,variableNames[i]));
             // getVariableObject searches both the data frame 'd' and the R environment
         if(variableObjects.back() != R_NilValue)
            variableTypes.push_back(getVariableType(variableObjects.back()));
         else {
            throw generalRbayzError("Variable name not found: "+variableNames[i]);
         }
      }
   }
   // Get the variance information, the part after V=, also first in a temporary string.
   // For a linear model specification (type 2) the whole part after '~' goes in varianceModel,
   // otherwise, the description must be a set of matrices and it is split on '*' and goes
   // in the vector of varianceNames.
   std::string tempvariance = getVarDescr(modelTerm);
   if (tempvariance != "") {
      if (tempvariance[0]=='~') {
         varianceType=3;
         varianceModel=tempvariance.substr(1,std::string::npos);
      }
      else {
         varianceNames = splitString(tempvariance,'*');
         varianceType=1;
         for(size_t i=0; i<varianceNames.size(); i++) {
            if( !(varianceNames[i]=="IDEN" || varianceNames[i]=="WEI" || varianceNames[i]=="MIXT") )
               varianceType=2; // if any variance-term is not IDEN, WEI, MIXT it is a correlated structure
         }
      }
   }
   // Get prior description
   priorModel = getPriorDescr(modelTerm);
}

