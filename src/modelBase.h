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
#include "nameTools.h"
#include "simpleMatrix.h"
#include "simpleVector.h"
#include "dcModelTerm.h"

class modelBase {
   
public:

   modelBase(dcModelTerm & modeldescr, modelBase * rmod);

   virtual ~modelBase() {
   }
   
   virtual void sample() = 0; 

   // prepForOutput is for model classes that need to make a transform of
   // parameters for output, the base class defines an 'empty' version.
   virtual void prepForOutput() { };

   std::string fname;
   simpleDblVector par, hpar;
   std::string parName, hparName;
   std::vector<std::string> parLabels, hparLabels, varNames;
   std::vector<Rcpp::RObject> varObjects;
   std::vector<int> varType;
   int hierType;
   std::string logPars; 

};

#endif /* modelBase_h */
