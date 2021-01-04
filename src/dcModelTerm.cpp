//  dcModelTerm.cpp

#include <Rcpp.h>
#include <vector>
#include <string>

dcModelTerm::dcModelTerm(std::string modelTerm, Rcpp::DataFrame &d) :
   varianceType(0), hierachType(0), varianceModel(""), priorModel(""), hierarchModel("")
{
   std::string tempnames = getVarNames(modelTerm);
   size_t pos = tempnames.find('/');      // '/' for hierarchical models
   if( (pos != std::string::npos) {
      size_t pos2 = tempnames.find('(');  // but proper models need parentheses somewhere
      if (pos2 == std::string::npos) {
         hierarchType=1;
         tempnames[pos]=':';
      }
      else {
         hierarchType=2;
      }
      std::string s1 = tempnames.substr(0,pos);
      std::string s2 =
      hasIndexVariable=TRUE;
   }
   varNames = splitString(tempnames,':');
   for(size_t i=0; i<varNames.size(); i++) {
      if (varNames[i]=="1" || varNames[i]=="0") {
         varObjects.push_back(R_NilValue);
         varType.push_back(0);
      }
      else {
         varObjects.push_back(getVariableObject(d,varNames[i]));
         if(varObjects.back() != R_NilValue)
            varType.push_back(getVariableType(varObjects.back()));
         else {
            throw generalRbayzError("Variable name not found: "+varNames[i]);
         }
      }
   }

}
varianceForm(0), varianceModel("")

   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
      int varianceType, hierarchType;
   std::vector<std::string> varianceNames;
   std::string varianceModel;
      std::string priorModel;
      std::string hierarchModel;
};

#endif /* dcModelTerm */
