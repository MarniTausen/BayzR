//
//  BayzR -- modelRanf.h
//  model class for random factor fit without correlation, derives from modelFactor.
//   - implements sample() method (see modelClasses.cpp for code)
//   - hpar is set as "var..." of size 1
//   - sample() also includes updating variance
//  Note: random effect with correation is quite different code, it derives from modelMatrix.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_h
#define modelRanf_h

#include <Rcpp.h>
#include "modelFactor.h"
#include "priorClasses.h"

class modelRanf : public modelFactor {

public:

   modelRanf(Rcpp::DataFrame &d, size_t col)  : modelFactor(d, col) {
      hpar.resize(1,1);
      hparName = "var." + parName;
   }

   ~modelRanf() {
   }

   void sample();

};

#endif /* modelRanf_h */
