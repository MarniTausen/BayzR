#include "modelBase.h"

// standard modelBase constructor
modelBase::modelBase(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
                 : par(), hpar(), gprior(modelTerm)
   {
      resid = e.data[2*resp];
      residPrec = e.data[2*resp+1];
      Nresid = e.nrow;
      parName = getVarNames(modelTerm);
      std::string tempnames = parName;
      size_t pos;
      // when making real hierarchical models, here the higher model could be isolated
      // and may be inserted as a new model inside the current model object?
      if( (pos=tempnames.find('/')) != std::string::npos) {
         hasIndexVariable=TRUE;
         tempnames[pos]=':';
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

