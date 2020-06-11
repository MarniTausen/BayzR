//
//  BayzR --- modelRanf_cor.hpp
//
// Computational class to to model ranf with correlation matrix (with V=K setting).
// This is a child-class of modelKernel because it works on a real matrix with regressions.
// - eigenvector data is already stored as coldata through constructor of modelTerm_realmat
// - eigenvalue data still needs to be set
// - hpar is one variance
// - par is set in modelTerm_realmat constructor and is size of coldata columns

//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_cor_h
#define modelRanf_cor_h

#include <Rcpp.h>
#include "modelKernel.h"

class modelRanf_cor : public modelKernel {

public:

   modelRanf_cor(Rcpp::DataFrame &d, size_t col) : modelKernel(d,col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3];
      hparName = "var." + parName;
   }

   ~modelRanf_cor() {
   }
   

   void sample() {
      dataKernel *K = static_cast<dataKernel *>(M);   // not yet a good solution for this,
                        // M is a pointer to dataMatrix but allocated to point to a dataKernel object.
                        // M can access the methods in dataMatrix, but not the extra eigenvalue
                        // information in dataKernel. Now K can get to the extra member variables.
      for(size_t k=0; k < K->nEvalUsed; k++) {
         resid_decorrect(k);
         collect_lhs_rhs(k);
         lhs = lhs + (1.0/( K->eval[k]*hpar[0] ));  // lhs with variance added
         par[k] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
         resid_correct(k);
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< K->nEvalUsed; k++)
         ssq += par[k]*par[k]/K->eval[k];
      hpar[0] = gprior.samplevar(ssq, K->nEvalUsed);
   }

};

#endif /* modelRanf_cor_h */
