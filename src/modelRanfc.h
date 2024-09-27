//  BayzR --- modelRanfc.h
//  Computational class to to model random effect with correlations (from using rn(..., V=K)).
//  Now also works with interactions and multiple kernels, but code now accepting ony max two kernels. 
//  Derives from modelFactor which sets-up the (interacting) factor data, the computational methods
//  though are more like modelMatrix and therefore parts of code are now same between modelMatrix and
//  ranfc. It would need multiple and virtual inheritance to re-use the modelMatrix methods here too.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanfc_h
#define modelRanfc_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "kernelMatrix.h"
#include "priorClasses.h"
#include "parsedModelTerm.h"

// update: ranfcor now also derive from modelFactor?
// the matrix data here was eigenvector data, it needs to come from a variance model now.

class modelRanfc : public modelFactor {

public:

   modelRanfc(parsedModelTerm & modeldescr, modelResp * rmod);
   ~modelRanfc() {}
   void sample();
   void accumFit(simpleDblVector & fit);
   void prepForOutput();
   std::vector<kernelMatrix*> kernelList;
   parVector *regcoeff;
   simpleDblVector fitval;
   std::vector<size_t> obsIndex;
   GenericPrior gprior;
   indepVarStr* varmodel;

};

#endif /* modelRanfc_h */
