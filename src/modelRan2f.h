//
//  BayzR --- modelRand2f
//  Computational class for two random factors.
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelRan2f_h
#define modelRan2f_h

#include <Rcpp.h>
#include "priorClasses.h"
#include "model2Factor.h"

class modelRan2f : public model2Factor {

public:

   modelRan2f(Rcpp::DataFrame &d, size_t col)  : model2Factor(d, col) {
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelRan2f() {
   }

   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=0; k<par.size(); k++) {
         par[k] = R::rnorm( (rhs[k]/(lhs[k]+(1/hpar[0]))), sqrt(1.0/(lhs[k]+(1/hpar[0]))));
      }
      resid_correct();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<par.size(); k++)
         ssq += par[k]*par[k];
      hpar[0] = gprior.samplevar(ssq,par.size());
   }

private:

   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.size(); k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs<intdata.size(); obs++) {
         k=intdata[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

protected:
   
   dataFactor *F1, *F2;

};

#endif /* modelRan2f_h */
