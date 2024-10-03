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
#include "parseFunctions.h"
#include "rbayzExceptions.h"

class modelRreg : public modelMatrix {

public:

   modelRreg(parsedModelTerm & modeldescr, modelResp * rmod)
         : modelMatrix(modeldescr, rmod)   {
      if(checkOptions(modeldescr.options, "V prior pvals")>0) {
         throw(generalRbayzError("ERROR: unrecognized option(s) in "+modeldescr.shortModelTerm));
      }
      if(modeldescr.options["pvals"]=="TRUE") {
         comp_frequentist_pvals = true;
         saveSamples = true;
         samplesFile = fopen();
      }
      // if pvals asked:
      // 1. may need to open a file to store sample info
      // 2. do we want a parVector to store pvalues? Can work well to get it in output, but this
      // requires a new model object, and if it is defined before chain is run, main() will
      // run sample(), parVector.collectStats() etc. on it.

   }

   ~modelRreg() {
   }
   
   void sample() {
      for(size_t k=0; k < M->ncol; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);   // update lhs and rhs variables
         lhs += varmodel->weights[k];
         par->val[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
      varmodel->sample();
      // need some thinking how to store sample info in file; likely every "saved" cycle.
      // main runs prepForOutput, and on a parVector main runs collectStats at the save intervals,
      // or it needs a new mechanism to switch on saving samples from output (which can be generic feature).
   }

   indepVarStr* varmodel;
   comp_frequentist_pvals=false;

};

// Here can start working on defining different variance structures for modelRreg
// For the moment, Rreg is only designed to accept variance models in the "indepVarStr" class, but this
// can include LASSO, BVS, DIAG / weighted, log-linear ...
class modelRregIden : public modelRreg {
public:
   modelRregIden(parsedModelTerm & pmdescr, modelResp * rmod)
      : modelRreg(pmdescr, rmod) {
      varmodel = new idenVarStr(pmdescr, this->par);
   }
};


#endif /* modelRreg */
