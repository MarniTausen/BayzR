//
//  rbayz -- model2Factor.h
//  Computational methods to work on two factors.
//
//  Created by Luc Janss on 29/05/2020.
//

#ifndef model2Factor_h
#define model2Factor_h

#include <Rcpp.h>
#include <cmath>
#include "modelBase.h"
#include "dataFactor.h"

class model2Factor : public modelBase {
   
public:
   
   model2Factor(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      F = new dataFactor(d, col);
      parLevelNames = F->labels;
      par.resize(parLevelNames.size(),0);
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~model2Factor() {
      delete F;
   }
   
protected:

   void resid_correct();
   void resid_decorrect();
   void collect_lhs_rhs();

   dataFactor *F1, *F2;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* model2Factor_h */
