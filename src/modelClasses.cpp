//  BayzR --- modelClasses.cpp
//

#include <Rcpp.h>
#include "modelRanfc.h"
#include "indexTools.h"

/* ------------------------- modelBase
*/

void modelBase::saveSamples() {

}

int modelBase::openSamplesFile() {

}

/* ------------------------- modelRanfc
*/

// old code that could only run cases with all known kernels up to max 2, and where for 2 kernels
// the kernels were merged in one combined kernel. The kernelList vector therefore remains size 1.
// It could be revived allowing for an option mergeKernels that merges known kernels in memory.
/*  old constructor
Ranfc_old(parsedModelTerm & modeldescr, modelResp * rmod)
{
   // For the moment all variance objects must be kernels
   for(size_t i=0; i<modeldescr.varianceObjects.size(); i++) {
      if (modeldescr.varianceObjects[i]==R_NilValue) {
         throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
      }
   }
   // Get the first kernel and then add (making kronecker products) with second etc., if available
   kernelList.push_back(new kernelMatrix(modeldescr.varianceObjects[0], modeldescr.varianceNames[0]));
   if (modeldescr.varianceNames.size()==2) {  // combine with a second kernel if present
      kernelMatrix* K2 = new kernelMatrix(modeldescr.varianceObjects[1], modeldescr.varianceNames[1]);
      kernelList[0]->addKernel(K2);
      delete K2;
   }
   if (modeldescr.varianceNames.size()>2) {  // need to think if I can keep combining kernels with addKernel()
      throw(generalRbayzError("Not yet ready to combine more than 2 kernels for interaction"));
   }
   // Here add a vector regcoeff (size K->ncol) to hold the regresssion on eigenvectors.
   // It is a parVector class so that the variance object can accept and work on it.
   regcoeff = new parVector(modeldescr, 0.0l, K->colnames);
   // obsIndex makes new level codes matching F->labels from every row in data to K->labels, it could
   // in principle replace the F->data and no need for the obsIndex vector.
   builObsIndex(obsIndex,F,K);
   fitval.initWith(F->nelem, 0.0l);
   // create the variance object - may need to move out as in ranfi
   varmodel = new idenVarStr(modeldescr, this->regcoeff);
}
*/

modelRanfc::modelRanfc(parsedModelTerm & modeldescr, modelResp * rmod)
           : modelFactor(modeldescr, rmod), regcoeff(), fitval(), gprior(modeldescr.options["prior"]) {
   // For the moment all variance objects must be kernels
   for(size_t i=0; i<modeldescr.varianceObjects.size(); i++) {
      if (modeldescr.varianceObjects[i]==R_NilValue) {
         throw(generalRbayzError("Mixing kernels with IDEN or other indep structures not yet possible"));
      }
   }
   for(size_t i=0; i<modeldescr.varianceObjects.size(); i++) {
      kernelList.push_back(new kernelMatrix(modeldescr.varianceObjects[i], modeldescr.varianceNames[i]));
   }
}

// sample method from old code
/*
Ranfc_old_sample() {
   // Update regressions on the eigenvectors
   double lhs, rhs;
   size_t matrixrow;
   double* colptr;
   kernelMatrix K = kernelList[0];
   for (size_t obs=0; obs < F->nelem; obs++)
      fitval[obs] = 0.0l;
   for(size_t col=0; col < K->ncol; col++) {
      colptr = K->data[col];
      // residual de-correction for this evec column
      for (size_t obs=0; obs < F->nelem; obs++)
         resid[obs] += regcoeff->val[col] * colptr[obsIndex[obs]];
      // Make the lhs and rhs and update this column regression
      lhs = 0.0l; rhs=0.0l;
      for (size_t obs=0; obs < F->nelem; obs++) {
         matrixrow = obsIndex[obs];
         rhs += colptr[matrixrow] * residPrec[obs] * (resid[obs]-fitval[obs]);
         lhs += colptr[matrixrow] * colptr[matrixrow] * residPrec[obs];
      }
      lhs += (1.0l / ( K->weights[col] * varmodel->par->val[0]) ) ;
      regcoeff->val[col] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
      // residual correction for this column with updated regression
      for (size_t obs=0; obs < F->nelem; obs++)
         fitval[obs] += regcoeff->val[col] * colptr[obsIndex[obs]];
   }
   for (size_t obs=0; obs < F->nelem; obs++)
      resid[obs] -= fitval[obs];
   // here need to replace with calling varmodel->sample(), but it needs
   // the diagVar structure.
   // update hyper-par (variance) using SSQ of random effects
   double ssq=0.0;
   for(size_t col=0; col< K->ncol; col++)
      ssq += regcoeff->val[col]*regcoeff->val[col]/K->weights[col];
   varmodel->par->val[0] = gprior.samplevar(ssq, K->ncol);
}
*/

// New sample() method with optimized scheme to compute statistics in blocks
void modelRanfc::sample() {
//   if(kernelList.size()==1) modelRanf::sample1();
//   else modelRanf::sample2();

}

void modelRanfc::accumFit(simpleDblVector & fit) {
   for (size_t obs=0; obs < F->nelem; obs++)
     fit[obs] += fitval[obs];
}

// prepForOutput puts the transform to random effects in the par-vector
// [!] this now only for 1, or 1 merged, kernel.
// Also: is this not the same as the fitted value??
void modelRanfc::prepForOutput() {
   kernelMatrix* K = kernelList[0];                  // kernelList is a std::vector<kernelMatrix*>
   for(size_t row=0; row< K->nrow; row++) {
      par->val[row]=0.0l;
      for(size_t col=0; col<K->ncol; col++) {
         par->val[row] += K->data[col][row] * regcoeff->val[col];
      }
   }
};
