//
//  rbayz -- modelFactor.h
//  Computational methods to work on one factor:
//  - declares and initialises a modelFactor object to work on
//  - sets sizes and names of parameter vectors - but not hpar, because that one differs
//    for derived classes
//  This is still not a concrete class -> see derived classes modelFixf and modelRanf.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelFactor_h
#define modelFactor_h

#include <Rcpp.h>
#include <cmath>
#include "modelBase.h"
#include "dataFactor.h"

class modelFactor : public modelBase {
   
public:
   
   modelFactor(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      F = new dataFactor(d, col);
      parLevelNames = F->labels;
      par.resize(parLevelNames.size(),0);
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~modelFactor() {
      delete F;
   }
   
protected:

   void resid_correct();
   void resid_decorrect();
   void collect_lhs_rhs();

   dataFactor *F;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* modelFactor_h */
