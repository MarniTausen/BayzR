//  dcModelTerm.cpp

#include <Rcpp.h>
#include <vector>
#include <string>
#include "dcModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

dcModelTerm::dcModelTerm(std::string modelTerm, Rcpp::DataFrame &d) :
   varianceType(0), varianceStruct(0), hierarchType(0), varianceLinMod(""), priorModel(""), hierarchModel("")
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
   // varianceType sets the main 3 cases of how variance is given:
   //   0: not given (there was no V= in the model term)
   //   1: linear model version (V=~....)
   //   2: list of structures (matrices) separated by stars e.g., V=VCOV[]*K1 or K1*K2*K3, etc
   // If varianceType is 1, the whole string after ~ is stored in VarianceLinMod.
   // If varianceType is 2, varianceNames, varianceParams and varianceObjects will have
   // one element for every component of the variance description (every part separated by stars) with
   //   varianceNames: the name such as "VCOV" or "K1"
   //   varianceParams: if there are parentheses after the name, the part that was inside
   //   varianceObjects: if the name can be found as an RObject, a link to the object
   // The name must be a known structure (now recognizing VCOV and MIXT) or name of an RObject,
   // otherwise an error is thrown.
   // Note: presence of not-Nil RObject means the variance name was not one of the bayz predefined
   // structures, and the RObject will be attempted to be interpreted as a kernel.
   // The varianceStruct sets if the variance structure can be handled by one of the INDEP
   // variance objects (varianceStruct=1) or if it needs a correlated structure (varianceStruct=2).
   // The linear model version, MIXT, weighted are INDEP, anything with VCOV of a kernel is CORR.
   // (Maybe need to expand for combinations VCOV*MIXT? - have to see how model objects are built...)
   std::string tempvariance = getVarDescr(modelTerm);
   if (tempvariance != "") {           // type 0 is set as default, only need to handle not zero
      if (tempvariance[0]=='~') {      // type 1
         varianceType=1;
         varianceStruct=1;
         varianceLinMod=tempvariance.substr(1,std::string::npos);
      }
      else {                           // must be (is interpreted as) type 2
         varianceType=2;
         std::vector<std::string> varianceElements = splitString(tempvariance,"*");
         for(size_t i=0; i<varianceElements.size(); i++) {
            size_t bracket = varianceElements[i].find_first_of("([");
            size_t len_tot = varianceElements[i].length();
            std::string name;
            std::string params;
            if (bracket == std::string::npos) {
               name = varianceElements[i];
               params = "";
            }
            else {
               size_t closeBrack = findClosingBrack(varianceElements[i], bracket);
               if(closeBrack != (len_tot-1) ) {
                  throw generalRbayzError("Unbalanced parentheses in: "+varianceElements[i]);
               }
               name = varianceElements[i].substr(0,bracket);
               params = varianceElements[i].substr(bracket+1,(len_tot-bracket-2));
            }
            varianceNames.push_back(name);
            varianceParams.push_back(params);
            if( name == "MIXT" ) {
               varianceObjects.push_back(R_NilValue);
               varianceStruct=1;
            }
            else {
               varianceObjects.push_back(getVariableObject(d,name));
               if(varianceObjects.back() == R_NilValue) {
                  throw generalRbayzError("Variance/kernel object not found in the R Environment: "
                       +varianceNames[i]);
               }
               varianceStruct=2;   // varianceStruct setting is not yet full-proof for multiple elements
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

