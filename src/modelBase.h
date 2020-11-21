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
      varNames = splitString(parName,':');
      // Search the varNames in the data-frame (column names) and fill varColIndex.
      // The searching of names is not efficient now, see comments in findDataColumn().
      for(size_t i=0; i<varNames.size(); i++) {
         if (varNames[i]=="1" || varNames[i]=="0")
            varColIndex.push_back(-1);
         else
            varColIndex.push_back(findDataColumn(d, varNames[i]));
      }
      // Also search the varNames in the R environment and fill varInEnvironment.
      // Already link it by setting an Robject? It would only be a reference, so
      // the overhead of having this Robject in every model object is minimal.
      // But .... needs to be a vector of Robjects?
      Rcpp::Environment Renv;
      for(size_t i=0; i<varNames.size(); i++) {
         if (varNames[i]=="1" || varNames[i]=="0")
            varInEnvironment.push_back(FALSE);
         else
            varInEnvironment.push_back(Renv.exists(varNames[i]));
      }
      // Check that every varName is found either in the data frame or in the environment
      // (Except for names "1" and "0"
      for(size_t i=0; i<varNames.size(); i++) {
         if( ! (varNames[i]=="1" || varNames[i]=="0") &&
                (varColIndex[i] == -1) && (varInEnvironment[i]==FALSE) )
            throw generalRbayzError("Variable name not found: "+varNames[i]);
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
   std::vector<int> varColIndex;
   std::vector<bool> varInEnvironment;
   
protected:
   double *resid, *residPrec;
   size_t Nresid;
   GenericPrior gprior;

};


#endif /* modelBase_h */
