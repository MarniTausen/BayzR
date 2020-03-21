//
//  BayzR -- modelFixf.h
//  model class for fixed factor fit, derives from modelFactor.
//   - implements sample() method (see modelClasses.cpp for code)
//   - hpar is not used in fixf model
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelFixf_h
#define modelFixf_h

#include <Rcpp.h>
#include <cmath>
#include "modelFactor.h"

class modelFixf : public modelFactor {
   
public:
   
   modelFixf(Rcpp::DataFrame &d, size_t col) : modelFactor(d, col) {
   }

   ~modelFixf() {
   }

   void sample();

};

#endif /* modelFixf_h */
