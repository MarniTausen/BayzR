// modelCoeff.h
// Parent of all coefficient models (all fix/ran/reg fitting objects).
// Created by Luc Janss on 08/03/2021.

#ifndef modelCoeff_h
#define modelCoeff_h

#include <Rcpp.h>
#include "modelResp.h"
#include "parsedModelTerm.h"

class modelCoeff : public modelBase {
   
public:

   modelCoeff(parsedModelTerm & modeldescr, modelResp * rmod) : modelBase() {
      // updating bayzR to handle hierarchical models can imply some (large?) changes here,
      // because then the "response model" can then be another coefficient model.
      // modelResp and modelCoeff are in different branches of the class hierarchy, and cannot
      // be cast from one to the other.
      // Maybe it needs a special interface class to help presenting a modelCoeff object
      // "as if" it is a modelResp class for the hierarchical model?
      Rcpp::Rcout << "In modelCoeff cstror\n";
      respModel=rmod;
      resid = respModel->resid.data;
      if(respModel->varModel==0) Rcpp::Rcout << "!respModel->varModel pointer is zero!\n";
      residPrec = respModel->varModel->weights.data;
      Nresid = respModel->par->nelem;
   }

   virtual ~modelCoeff() {
   }

   virtual void accumFit(simpleDblVector & fit) = 0;

   modelResp* respModel;
   // The following is for convenience so that all modelCoeff objects have direct
   // pointers to residuals and residual variance, and it only needs to be set once
   // here in modelCoeff cstr'or.
   double *resid=NULL, *residPrec=NULL;
   size_t Nresid=0;

};

#endif /* modelCoeff_h */
