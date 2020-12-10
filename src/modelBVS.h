//
//  BayzR --- modelBVS.h
//  Bayesian Variable Selection model - derives from modelMatrix to be a 'sibling' model
//  to modelRreg and work on the same data structures as Rreg.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelBVS_h
#define modelBVS_h

#include <Rcpp.h>
#include "modelMatrix.h"

class modelBVS : public modelMatrix {

public:

   modelBVS(std::string modelTerm, Rcpp::DataFrame &d, simpleMatrix &e, size_t resp)
         : modelRreg(modelTerm, d, e, resp)
   {
      par = regcoeff;    // same as in Rreg
      hpar.initWith(par.nelem+2,1.0l);
      hparName[0] = "pi0." + parName;
      hparName[1] = "tau0." + parName;
   }

   ~modelBVS() {
   }
   
   void sample() {
      update_regressions();
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]/M->weights[k];
      hpar[0] = gprior.samplevar(ssq, M->ncol);
   }


};

#endif /* modelBVS */
