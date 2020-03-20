//
//  modelTerm_ran_cor.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelTerm_ran_cor_h
#define modelTerm_ran_cor_h

#include <Rcpp.h>
#include "modelTerm_realmat.h"
#include "priorClasses.h"
#include "parseColNames.h"
#include "data_kernel.h"

// random-correlated model term, built from ranf() model term with V=K setting.
// This is a child-class of modelTerm_realmat because it works on a real matrix with regressions.
// - eigenvector data is already stored as coldata through constructor of modelTerm_realmat
// - eigenvalue data still needs to be set
// - hpar is one variance
// - par is set in modelTerm_realmat constructor and is size of coldata columns

class modelTerm_ran_cor : public modelTerm_factor {

public:

   modelTerm_ran_cor(Rcpp::DataFrame &d, size_t col) : v((Rcpp::Robject)(d[col])) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3];
      hparName = "var." + parName;
   }

   ~modelTerm_ran_cor() {
   }
   

   void sample() {
      for(size_t k=0; k<nEvalUsed; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         lhs = lhs + (1.0/(eval[k]*hpar[0]));  // lhs with variance added
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<nEvalUsed; k++)
         ssq += par[k]*par[k]/eval[k];
      hpar[0] = gprior.samplevar(ssq,nEvalUsed);
   }
   
protected:
   // These methods may be general for working on matrices, it could be organised
   // by making a class that handles factor+matrix computation, or with multiple inheritance.
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

private:

   dataKernel v;
   Rcpp::NumericVector eval;
   Rcpp::IntegerVector update;
   
};

#endif /* modelTerm_ran_cor_h */
