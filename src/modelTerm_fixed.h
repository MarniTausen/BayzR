//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_fixed_h
#define modelTerm_fixed_h

#include <Rcpp.h>
#include <cmath>
#include "modelTerm_factor.h"

// Model-term for fixed effect (an R factor)
// - the coldata is IntegerVector
// - par vector is size number of levels of the factor
// - hpar is not used here

class modelTerm_fixed : public modelTerm_factor {
   
public:
   
   modelTerm_fixed(Rcpp::DataFrame &d, size_t col) : modelTerm_factor(d, col) {
   }
   
   ~modelTerm_fixed() {
   }
   
   void sample() {
      resid_decorrect();
      collect_lhs_rhs();
      for(size_t k=1; k<par.size(); k++) {  // par[0] remains zero!, this runs from k=1
         if (lhs[k]>0)
            par[k] = R::rnorm( (rhs[k]/lhs[k]), sqrt(1.0/lhs[k]));
         else
            par[k]=0.0;
      }
      resid_correct();
   }

private:

};

#endif /* modelTerm_fixed_h */
