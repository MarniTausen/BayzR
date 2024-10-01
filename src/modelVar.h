// Base class for variance models
// Defines that that all modelVar objects have a gprior and parVector* coeff
// (points to par from coefficient model) to work on.
//
// Created by Luc Janss on 03/03/2021.

#ifndef modelVar_h
#define modelVar_h

#include <Rcpp.h>
#include <vector>
#include "modelBase.h"
#include "priorClasses.h"
#include "parsedModelTerm.h"
//#include <unistd.h>

class modelVar : public modelBase {
   
public:

   modelVar(parsedModelTerm & modeldescr)
      : modelBase(), gprior(modeldescr.options["prior"]) {
   }

   virtual ~modelVar() { }
   
protected:
   GenericPrior gprior;
   parVector* coefpar=0;

};

#endif /* modelVar_h */