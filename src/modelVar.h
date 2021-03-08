// Base class for variance models
// Created by Luc Janss on 03/03/2021.

#ifndef modelVar_h
#define modelVar_h

#include <Rcpp.h>
#include <vector>
#include "modelBase.h"
#include "priorClasses.h"
#include "simpleVector.h"
#include "dcModelTerm.h"

class modelVar : public modelBase {
   
public:

   modelVar(dcModelTerm & modeldescr, modelBase * rmod)
      : modelBase(modeldescr, rmod), gprior(modeldescr.priorModel) { }

   virtual ~modelVar() { }
   
//   simpleDblVector par, hpar;
//   std::string parName, hparName;
//   std::vector<std::string> parLabels, hparLabels;

protected:
   GenericPrior gprior;

};

#endif /* modelVar_h */