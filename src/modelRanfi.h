// BayzR -- modelRanfi.h
// Model class for RANdom Factor with Indep variance structures.
// Derives from modelFactor and has pointer to IndepVarStr object to model the variance.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanfi_h
#define modelRanfi_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "indepVarStr.h"
//#include <unistd.h>

class modelRanfi : public modelFactor {

public:

   modelRanfi(parsedModelTerm & modeldescr, modelResp * rmod)
      : modelFactor(modeldescr, rmod) {
   }

   ~modelRanfi() {
      delete varmodel;
   }

   void sample() {
      // sample() method for random effects with indep var-structure: the data
      // corrections and LHS and RHS can be made using the methods from modelFactor.
      // To update parameters as random effects varmodel->weights[k] are added in LHS.
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=0; k<par->nelem; k++) {
         lhs[k] += varmodel->weights[k];
         par->val[k] = R::rnorm((rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
      }
      resid_correct();
   }

   void sampleHpars() {
      varmodel->sample();
   }

   indepVarStr* varmodel;

};

class modelRanFacIden : public modelRanfi {
public:
   modelRanFacIden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRanfi(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};


#endif /* modelRanfi_h */
