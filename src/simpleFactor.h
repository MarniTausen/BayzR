//  R/Bayz
//  simpleFactor.h - there is no 'native' factor in Rcpp, so I created one for use inside
//     the C++ code using simpleIntVector with added labels (the R 'levels').
//     In this shape it is not exportable/useable as a proper R factor, but it could be
//     quite easy to transform if needed.
//
//  Created by Luc Janss on 11/03/2024.
//

#ifndef simpleFactor_h
#define simpleFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "simpleVector.h"

class simpleFactor : public simpleIntVector {
public:
   simpleFactor(Rcpp::RObject col, std::string name);
   ~simpleFactor() {
   } ;
   std::vector<std::string> labels;
   std::string name;
   std::vector<std::string> back2vecstring();
};

#endif /* simpleFactor_h */
