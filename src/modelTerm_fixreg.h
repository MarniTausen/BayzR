//
//  modelTerm_fixreg.hpp
//  rbayz
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelTerm_fixreg_h
#define modelTerm_fixreg_h

#include <Rcpp.h>
#include <cmath>
#include "modelTerm_realvec.h"

// Model-term for fixed regression (built from freg() in the model).
// This class derives from modelTerm_realvec, which
// has defined storaged and methods to collect statistics for a real data vector.
// - the coldata is NumericVector
// - par vector is size 1, we only update par[0] here
// - hpar is not used here

class modelTerm_fixreg : public modelTerm_realvec {
   
public:
   
   modelTerm_fixreg(Rcpp::DataFrame &d, size_t col) : modelTerm_realvec(d, col) {
   }
   
   ~modelTerm_fixreg() {
   }
   
   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      par[0] = R::rnorm( (rhs/lhs), sqrt(1.0/lhs));
      resid_correct();
   }

private:

};

#endif /* modelTerm_fixreg_h */
