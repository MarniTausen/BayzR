//
//  BayzR -- modelFreg.h
//  Computational class defining methods to work on one covariate.
//  - the parameter vector has length 1, and no levelNames.
//  - hpar is not used.
//  - also defines sample(). So far only one model-object that uses one
//    covariate, so it is combining all de/correction, lhs/rhs statistics
//    and sample code.
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelFreg_h
#define modelFreg_h

#include <Rcpp.h>
#include <cmath>
#include "modelBase.h"
#include "dataCovar.h"

class modelFreg : public modelBase {
   
public:
   
   modelFreg(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      C = new dataCovar(d, col);
      par.resize(1,0);
   }
   
   ~modelFreg() {
      delete C;
   }
   
   void sample();
   
protected:

   void resid_correct();
   void resid_decorrect();
   void collect_lhs_rhs();

   dataCovar *C;
   double lhs, rhs;
   std::vector<double> fit;  // not used yet
   
};

#endif /* modelFreg_h */
