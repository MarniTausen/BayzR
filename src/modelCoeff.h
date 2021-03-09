// modelCoeff.h
// Parent of all coefficient models (all fix/ran/reg fitting objects).
// Created by Luc Janss on 08/03/2021.

#ifndef modelCoeff_h
#define modelCoeff_h

#include <Rcpp.h>
#include "modelResp.h"
#include "dcModelTerm.h"

class modelCoeff : public modelBase {
   
public:

   modelCoeff(dcModelTerm & modeldescr, modelBase * rmod) : modelBase(modeldescr, rmod) {
      respModel = dynamic_cast<modelResp *>(rmod);  // for now rmod passed is a modelResp, but maybe not
                                                    // when using hierarchies -> rmod can become a modelCoeff
      if (modeldescr.hierarchType == 0 || modeldescr.hierarchType == 1) {
         if (respModel != NULL) {
            resid = respModel->par.data;
            residPrec = respModel->varModel->weights.data;
            Nresid = respModel->par.nelem;
         }
      } else { // Hierarchical model:
               // don't know yet how to do here, there can be different cases,
               // but the rmod may not have a resid and residPrec vector, and the rmod
               // will not be a respModel object ....
//         respModel = rmod;
      }
   }

   virtual ~modelCoeff() {
   }
   
   modelResp* respModel;
   double *resid=NULL, *residPrec=NULL;
   size_t Nresid=0;

};

#endif /* modelCoeff_h */
