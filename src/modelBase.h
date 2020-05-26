//
//  modelBase.h
//
//  Base class for model (computational) classes.
//   - declares and set the resid and residPrec vector to the one-by-last and last
//     column of the model-frame
//   - declare par and hpar vectors, and vectors to hold names of the parameter and
//     hyper-parameter (parName, hparName), and labels for their levels ((h)parLevelNames).
//    * the base class cannot set sizes or fill all names, if a derived class fails to
//       intialize properly problems can occur!
//    * for parName a default is set, but sometimes needs changing in a derived class.
//   - gprior member variable is declared and filled from information in the data column
//   - define sample() as pure virtual function, all derived classes should make an
//     implementation of this for proper working!
//
//  Created by Luc Janss on 03/08/2018.

#ifndef modelBase_h
#define modelBase_h

#include <Rcpp.h>
#include <vector>
#include "parseColNames.h"
#include "priorClasses.h"

class modelBase {
   
public:

   modelBase(Rcpp::DataFrame &d, size_t col) : gprior(d, col) {
      resid = d[d.size()-2];
      residPrec = d[d.size()-1];
      if(getWrapName(d,col)=="")
         parName = getVarName(d, col);
      else
         parName = getWrapName(d, col) + "." + getVarName(d, col);
   }

   virtual ~modelBase() {
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


#endif /* modelBase_h */
