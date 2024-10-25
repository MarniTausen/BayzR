//
//  BayzR -- modelFreg.h
//  Computational class defining methods to work on one covariate.
//  - the parameter vector has length 1, and no levelNames.
//  - no hyper parameters - sampleHpar is empty.
//  - also defines sample(). So far only one model-object that uses one
//    covariate, so it is combining all de/correction, lhs/rhs statistics
//    and sample code.
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelFreg_h
#define modelFreg_h

#include <Rcpp.h>
#include <cmath>
#include "modelCoeff.h"
#include "dataCovar.h"
#include "parsedModelTerm.h"

class modelFreg : public modelCoeff {
   
public:
   
   modelFreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelCoeff(modeldescr, rmod)
   {
      if(modeldescr.variableTypes[0] != 3)
         throw generalRbayzError("Not a numeric vector input: "+modeldescr.variableNames[0]);
      C = new dataCovar(modeldescr.variableObjects[0]);
      par = new parVector(modeldescr, 0.0l);
      par->traced=1;
   }
   
   ~modelFreg() {
      delete C;
   }
   
   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      par->val[0] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
      resid_correct();
   }

   void sampleHpars() {}

   void accumFit(simpleDblVector & fit) {
      for (size_t obs=0; obs < C->nelem; obs++)
        fit[obs] += par->val[0] * C->data[obs];
   }

   void resid_correct() {
      for (size_t obs=0; obs < C->nelem; obs++)
        resid[obs] -= par->val[0] * C->data[obs];
   }

   void resid_decorrect() {
      for (size_t obs=0; obs < C->nelem; obs++)
        resid[obs] += par->val[0] * C->data[obs];
   }

   void collect_lhs_rhs() {
      lhs = 0.0; rhs=0.0;
      for (size_t obs=0; obs < C->nelem; obs++) {
         rhs += residPrec[obs] * resid[obs] * C->data[obs];
         lhs += C->data[obs] * residPrec[obs] * C->data[obs];
      }
   }

   dataCovar *C;
   double lhs, rhs;
   std::vector<double> fit;  // not used yet
   
};

#endif /* modelFreg_h */
