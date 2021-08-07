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
   size_t pos1, pos2, pos3;
   if ( (pos1 = modelTerm.find('(')) == std::string::npos)
      funcName = "";
   else
      funcName = modelTerm.substr(0, pos1);
   // Variables name(s) need some pre-treatment to handle hierarchy specifications
   std::string tempnames = getVarNames(modelTerm);
   pos1 = tempnames.find('/');
   pos2 = tempnames.find('(');
   pos3 = tempnames.find('|');
   if( pos1 != std::string::npos && pos2 != std::string::npos ) {
      // both '/' and '(' : hierarchical model. Remove the hierarchical part and
      // put it in hierarchModel string, all processing continuous without it.
      // The string allVariableNames will also not have the hierarchy part.
      hierarchType=2;
      hierarchModel = tempnames.substr(pos1+1, std::string::npos);
      tempnames.erase(pos1, std::string::npos);
   }
   // Here insert finding the variable patterns based on :|/
   // ....
   allVariableNames = tempnames;
   variableNames = splitString(tempnames,":|/");
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
            throw generalRbayzError("Variable object not found: "+variableNames[i]);
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
         // this may need reviewing or extension later:
         // 1) not sure how to handle or to accept variance with multiple indep structures, e.g. WEI*MIXT etc.
         // 2) WEI and maybe also MIXT may get extra parameters or variables
         varianceNames = splitString(tempvariance,"*");
         varianceType=1;
         for(size_t i=0; i<varianceNames.size(); i++) {
            if( varianceNames[i]=="IDEN" || varianceNames[i]=="WEI" || varianceNames[i]=="MIXT" ) {
               varianceObjects.push_back(R_NilValue);
            }
            else {
               varianceObjects.push_back(getVariableObject(d,varianceNames[i]));
               if(varianceObjects.back() == R_NilValue) {
                  throw generalRbayzError("Kernel object not found: "+varianceNames[i]);
               }
               varianceType=2; // flag variance type 2 if some matrices are not IDEN, WEI, MIXT
            }
         }
      }
   }
   // Get prior description
   priorModel = getPriorDescr(modelTerm);
   // Check / get a logging option
   std::string s = getOptionText(modelTerm, "log=");
   if (s == "")    // if s is empty, set logging to "def" (default)
      logging = "def";
   else
      logging=s;
   // can add checks that the logging setting is from an allowed set of settings
}

