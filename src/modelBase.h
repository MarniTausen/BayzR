//
//  modelBase.h
//
//  Base class for model (computational) classes.
//  Of course, this constructor is run first for construction of every model-class object,
//  so it allows to prepare things that can be used by all model-classes.
//  - resid and residPred pointers are set, they depend on the response that is handled
//  - define all common variables: parameters, names, parameter labels
//  - this also checks and finds/links variable-names to data-frame columns, or
//    checks that they can be found in the R environment
//  - defines that all model objects have a grprior objects and retrieves information
//    from a prior= setting in it, when available.
//  - defines sample() as pure virtual function, all derived classes should make an
//     implementation of this for proper working!
// Note:
//    * the base class cannot set sizes or fill all names of paranmeters, if a derived
//      class fails to intialize properly problems can occur!
//    * for parName a default is set as the variable name, but sometimes needs changing.
//
//  Created by Luc Janss on 03/08/2018.

#ifndef modelBase_h
#define modelBase_h

#include <Rcpp.h>
#include <vector>
#include "parseFunctions.h"
#include "priorClasses.h"
#include "nameTools.h"
#include "simpleMatrix.h"

class modelBase {
   
public:

   modelBase(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
                 : gprior(modelTerm)
   {
      resid = e.data[2*resp];
      residPrec = e.data[2*resp+1];
      Nresid = e.nrow;
      parName = getVarNames(modelTerm);
      size_t pos;
      // when making real hierarchical models, here the higher model could be isolated
      // and may be inserted as a new model inside the current model object?
      if( (pos=parName.find(';')) != std::string::npos) {
         hasIndexVariable=TRUE;
         parName[pos]=':';
      }
      varNames = splitString(parName,':');
      for(size_t i=0; i<varNames.size(); i++) {
         if (varNames[i]=="1" || varNames[i]=="0") {
            varObjects.push_back(R_NilValue);
            varType.push_back(0);
         }
         else {
            varObjects.push_back(getVariableObject(d,varNames[i]);
            if(varObjects.back() != R_NilValue)
               varType.push_back(getVariableType(varObjects.back()));
            else {
               throw generalRbayzError("Variable name not found: "+varNames[i]);
            }
         }
      }
   }

   virtual ~modelBase() {
   }
   
   // sample is pure virtual, the base class cannot give a definition, and
   // all derived classes must implement it to create a concrete class.
   virtual void sample()=0;

   // prepForOutput is for model classes that need to make a transform of
   // parameters for output, the base class defines an 'empty' version.
   virtual void prepForOutput() { };
   
   std::vector<double> par, hpar;
   std::string parName, hparName;
   std::vector<std::string> parLabels, hparLabels, varNames;
   std::vector<Rcpp::RObject> varObjects;
   std::vector<int> varType;
   bool hasIndexVariable=FALSE;
   
protected:
   double *resid, *residPrec;
   size_t Nresid;
   GenericPrior gprior;

};

#endif /* modelBase_h */
