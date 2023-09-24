//
//  BayzR --- modelRreg.hpp
//
//  Computational class to model random regressions on a matrix of covariates.
//  This uses methods from modelMatrix, only some parameter names need to be set.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRreg_h
#define modelRreg_h

#include <Rcpp.h>
#include "modelMatrix.h"
#include "indepVarStr.h"
#include "dataMatrix.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelMatrix(modeldescr, rmod)   {
   }

   ~modelRreg() {
   }
   
   void sample() {
      for(size_t k=0; k < M->ncol; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);   // update lhs and rhs variables
         lhs += varmodel->weights[k];
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
   }

   indepVarStr* varmodel;

};

// Here can start working on defining different variance structures for modelRreg

class modelRregIden : public modelRreg {
public:
   modelRregIden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};


#endif /* modelRreg */
