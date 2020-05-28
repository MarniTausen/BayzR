//
//  BayzR -- modelClasses.cpp
//
//  Created by Luc Janss on 24/01/2020.
//

#include <Rcpp.h>
#include <cmath>
#include "modelResp.h"
#include "modelMean.h"
#include "modelLiab.h"
#include "modelFactor.h"
#include "modelFixf.h"
#include "modelRanf.h"
#include "modelFreg.h"
#include "modelMatrix.h"

void modelResp::sample() {
   // Continuous data: sample() is only updating residual variance
   size_t obs, nobs=resid.size();
   double sum=0.0;
   for (obs=0; obs<nobs; obs++)
      sum += resid[obs]*resid[obs]*residPrec[obs];
   sum *= hpar[0];  // the sum was computed divided by the old variance!
   hpar[0] = gprior.samplevar(sum, nobs);
   // the residPrec vector should be re-filled
   double temp = 1.0/hpar[0];
   for (obs=0; obs<nobs; obs++)
      residPrec[obs] = temp;
}

void modelMean::sample() {
   size_t obs, nobs=resid.size();
   double sum=0.0, temp=0.0;
   for (obs=0; obs<nobs; obs++) resid[obs] += par[0];
   for (obs=0; obs<nobs; obs++) {
      sum += resid[obs]*residPrec[obs];
      temp += residPrec[obs];
   }
   par[0] = R::rnorm((sum/temp), sqrt(1.0/temp));
   for (obs=0; obs<nobs; obs++) resid[obs] -= par[0];
}

std::pair<double, double> modelLiab::liabBounds(int i) {
   std::pair<double, double> bounds;
   // not yet finished here ...
   return bounds;
}

void modelFactor::resid_correct() {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] -= par[F->data[obs]];
}

void modelFactor::resid_decorrect() {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] += par[F->data[obs]];
}

void modelFactor::collect_lhs_rhs() {
   size_t k;
   for(k=0; k<par.size(); k++) {
      rhs[k] = 0.0;
      lhs[k] = 0.0;
   }
   for (size_t obs=0; obs < F->data.size(); obs++) {
      k=F->data[obs];
      rhs[k] += residPrec[obs] * resid[obs];
      lhs[k] += residPrec[obs];
   }
}

void modelFixf::sample() {
   resid_decorrect();
   collect_lhs_rhs();
   for(size_t k=1; k<par.size(); k++) {  // in fixf par[0] remains zero!, this runs from k=1
      if (lhs[k]>0)                      // if lhs (inverse variance) is zero estimate will be set to 0
         par[k] = R::rnorm( (rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
      else
         par[k]=0.0;
   }
   resid_correct();
}

void modelRanf::sample() {
   resid_decorrect();
   collect_lhs_rhs();
   for(size_t k=0; k<par.size(); k++) {
      // random effect: add 1/hpar[0] in lhs
      par[k] = R::rnorm( (rhs[k]/(lhs[k]+(1/hpar[0]))), sqrt(1.0/(lhs[k]+(1/hpar[0]))));
   }
   resid_correct();
   // update hyper-par (variance) using SSQ of random effects
   double ssq=0.0;
   for(size_t k=0; k<par.size(); k++)
      ssq += par[k]*par[k];
   hpar[0] = gprior.samplevar(ssq,par.size());
}

void modelFreg::resid_correct() {
   for (size_t obs=0; obs < C->data.size(); obs++)
      resid[obs] -= par[0] * C->data[obs];
}

void modelFreg::resid_decorrect() {
   for (size_t obs=0; obs < C->data.size(); obs++)
      resid[obs] += par[0] * C->data[obs];
}

void modelFreg::collect_lhs_rhs() {
   lhs = 0.0; rhs=0.0;
   for (size_t obs=0; obs < C->data.size(); obs++) {
      rhs += residPrec[obs] * resid[obs] * C->data[obs];
      lhs += C->data[obs] * residPrec[obs] * C->data[obs];
   }
}

void modelFreg::sample() {
   resid_decorrect();
   collect_lhs_rhs();
   par[0] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
   resid_correct();
}

void modelMatrix::resid_correct(size_t col) {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] -= par[col] * M->data(F->data(obs),col);
}

void modelMatrix::resid_decorrect(size_t col) {
   for (size_t obs=0; obs < F->data.size(); obs++)
      resid[obs] += par[col] * M->data(F->data(obs),col);
}

void modelMatrix::collect_lhs_rhs(size_t col) {
   lhs = 0.0; rhs=0.0;
   size_t rowlevel;
   for (size_t obs=0; obs < F->data.size(); obs++) {
      rowlevel = F->data(obs);
      rhs += M->data(rowlevel,col) * residPrec[obs] * resid[obs];
      lhs += M->data(rowlevel,col) * M->data(rowlevel,col) * residPrec[obs];
   }
}

/* old code from ranf_cor for comparison
 void resid_correct(size_t col) {
    for (size_t obs=0; obs < coldata.size(); obs++)
       resid[obs] -= par[col] * matrixdata(coldata(obs),col);
 }
 
 void resid_decorrect(size_t col) {
    for (size_t obs=0; obs < coldata.size(); obs++)
       resid[obs] += par[col] * matrixdata(coldata(obs),col);
 }

 void collect_lhs_rhs(size_t col) {
    lhs = 0.0; rhs=0.0;
    size_t rowlevel;
    for (size_t obs=0; obs < coldata.size(); obs++) {
       rowlevel = coldata(obs);
       rhs += matrixdata(rowlevel,col) * residPrec[obs] * resid[obs];
       lhs += matrixdata(rowlevel,col) * matrixdata(rowlevel,col) * residPrec[obs];
    }
 }

 */
