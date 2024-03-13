//
//  rbayz -- modelFactor.h
//  Computational methods to work on one factor that is modelled fixed, or random without
//  correlations (modelFixf and modelRanf derive from this and only use small modifications
//  to use common code from modelFactor).
//  - declares and initialises a modelFactor object
//  - sets sizes and names of parameter vectors - but not hpar, because that one differs
//    for derived classes (fixf has no hpar)
//  - now sets up factor with any number of interactions using new features from dataFactor
//  This is still not a concrete class -> see derived classes modelFixf and modelRanf.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelFactor_h
#define modelFactor_h

#include <Rcpp.h>
#include <cmath>
#include "modelCoeff.h"
#include "dataFactor.h"

class modelFactor : public modelCoeff {
   
public:
   
   modelFactor(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod)
   {
      F = new dataFactor(modeldescr.variableObjects, modeldescr.variableNames);
      par = new parVector(modeldescr, 0.0l, F->labels);
      lhs.resize(F->labels.size(),0);
      rhs.resize(F->labels.size(),0);
   }
   
   ~modelFactor() {
      delete F;
   }

   void accumFit(simpleDblVector & fit) {
      for (size_t obs=0; obs < F->nelem; obs++)
        fit[obs] += par->val[F->levcode[obs]];
   }

   
protected:

   void resid_correct() {
      for (size_t obs=0; obs < F->nelem; obs++)
        resid[obs] -= par->val[F->levcode[obs]];
   }

   void resid_decorrect() {
      for (size_t obs=0; obs < F->nelem; obs++)
        resid[obs] += par->val[F->levcode[obs]];
   }

   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par->nelem; k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs < F->nelem; obs++) {
         k=F->levcode[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

   dataFactor *F;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations
                                          // maybe faster using the simpleVector class?

};

#endif /* modelFactor_h */
