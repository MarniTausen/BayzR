//
//  rbayz -- modelFactor.h
//  Computational methods to work on one factor.
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

   void resid_correct() {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] -= par[F->data[obs]];
   }
   
   void resid_decorrect() {
      for (size_t obs=0; obs < F->data.size(); obs++)
         resid[obs] += par[F->data[obs]];
   }

   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.size(); k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs < F->data.size(); obs++) {
         k=F->data[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

   dataFactor *F;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* modelFactor_h */
