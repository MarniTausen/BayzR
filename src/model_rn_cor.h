//  BayzR --- model_rn_cor.h
//  Computational classes to model random effect with correlations (from using rn(..., V=K)).
//  Now also works with interactions and multiple kernels, but code now accepting ony max two kernels. 
//  Derives from modelFactor which sets-up the (interacting) factor data, the computational methods
//  though are more like modelMatrix and therefore parts of code are now same between modelMatrix and
//  ranfc. It would need multiple and virtual inheritance to re-use the modelMatrix methods here too.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef model_rn_cor_h
#define model_rn_cor_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "kernelMatrix.h"
#include "priorClasses.h"
#include "parsedModelTerm.h"

// update: ranfcor now also derive from modelFactor?
// the matrix data here was eigenvector data, it needs to come from a variance model now.

class model_rn_cor_k0 : public modelFactor {

public:

   model_rn_cor_k0(parsedModelTerm & modeldescr, modelResp * rmod);
   ~model_rn_cor_k0() {}
   void sample();
   void sampleHpars();
   void accumFit(simpleDblVector & fit);
   void prepForOutput();
   std::vector<kernelMatrix*> kernelList;
   parVector *regcoeff;
   simpleDblVector fitval;
   std::vector<size_t> obsIndex;
   indepVarStr* varmodel;

};



#endif /* model_rn_cor_h */
