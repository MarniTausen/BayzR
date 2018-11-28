//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_h
#define modelTerm_h

#include <Rcpp.h>
#include <vector>
#include "parseColNames.h"
#include "priorClasses.h"

/********
 Base class for model-terms.
 
 All model-terms derive from the same base-class 'modelTerm'.
 This allows in the main function to use vector<modelTerm *> to collect pointers to all model-terms, and
 do things like calling all sample() methods on all modelTerm objects.
 
 The base-class can organize some common things for all model-terms:
 - declare and set the resid and residPrec vector, to the one-by-last and last column of the model-frame
 - declare and set the parName to default column-WrapName + "." + column-VarName, however, some model-terms may
   need to redo/modify it when it needs to be something else. There is no default hparName set in the base-class.
 - declare object of GenericPrior class and retrieve possible prior information in it (may be empty)
 
 Each derived class should do:
 - declare datacol (types differ between model-terms) and point it to the column in the model-frame that
   the object works on
 - set sizes used of par and hpar vector (this can remain zero if not used)
 - if needed, reset parName, and set hparName
 - construct / retrieve parLevelNames and hparLevelNames
 - implement sample() method
 
 Some special parts:
 All model-terms work on / use the common resid and residPrec vector, which are added as two extra columns
 in the model-frame. The modelTerm_resp object fills the initial resid and residPrec vectors in its constructor,
 and also updates these vectors in its sample() method.
 
********/

class modelTerm {
   
public:

   modelTerm(Rcpp::DataFrame &d, size_t col) : gprior(d, col) {
      resid = d[d.size()-2];
      residPrec = d[d.size()-1];
      parName = getWrapName(d, col) + "." + getVarName(d, col);
   }

   virtual ~modelTerm() {
   }
   
   virtual void sample()=0;
   
   std::vector<double> par, hpar;
   std::string parName, hparName;
   Rcpp::CharacterVector parLevelNames, hparLevelNames;
   
protected:
   Rcpp::NumericVector resid;
   Rcpp::NumericVector residPrec;
   GenericPrior gprior;

};


#endif /* modelTerm_h */
