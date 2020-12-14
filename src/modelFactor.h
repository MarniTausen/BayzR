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
#include "modelBase.h"
#include "dataFactor.h"

void CharVec2cpp(std::vector<std::string> & labels, Rcpp::CharacterVector templabels);

class modelFactor : public modelBase {
   
public:
   
   modelFactor(std::string modelTerm, Rcpp::DataFrame &d, modelBase * rmod)
         : modelBase(modelTerm, d, rmod)
   {
      F = new dataFactor();
      for(size_t i=0; i<varNames.size(); i++) {
         if (varType[i] != 1)
            throw generalRbayzError("Variable is not a factor: "+varNames[i]);
         F->addVariable(varObjects[i]);
      }

      parLabels = F->labels;            // I think this would make a copy, but actually only
                                        // a reference is enough. The double storage arises here
                                        // because all model-classes have a parLabels vector, which
                                        // here happens to be equal to the labels of the Factor object.
      par.initWith(parLabels.size(),0.0l);
      lhs.resize(parLabels.size(),0);
      rhs.resize(parLabels.size(),0);
   }
   
   ~modelFactor() {
      delete F;
   }
   
protected:

   void resid_correct() {
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] -= par[F->data[obs]];
   }

   void resid_decorrect() {
      for (size_t obs=0; obs < F->data.nelem; obs++)
         resid[obs] += par[F->data[obs]];
   }

   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.nelem; k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs < F->data.nelem; obs++) {
         k=F->data[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

   dataFactor *F;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* modelFactor_h */
