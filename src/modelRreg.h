//
//  BayzR --- modelRreg.hpp
//
//  Computational class to to model random regressions on a matrix of covariates.
//  This used all methods from modelMatrix, only some parameter names need to be set.
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

   modelRreg(dcModelTerm & modeldescr, modelBase * rmod)
         : modelMatrix(modeldescr, rmod)
   {
      if (M->colnames.size() >0)    // if there are no colnames in matrix, default ones
         parLabels = M->colnames;   // will be inserted in the parameter names list
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

#endif /* modelRreg */
