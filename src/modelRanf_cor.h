//
//  BayzR --- modelRanf_cor.hpp
//
//  Computational class to to model ranf with correlation matrix (with V=K setting).
//  This used all methods from modelMatrix, only some parameter names need to be set.
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelRanf_cor_h
#define modelRanf_cor_h

#include <Rcpp.h>
#include "modelMatrix.h"
#include "dataMatrix.h"

class modelRanf_cor : public modelMatrix {

public:

   modelRanf_cor(Rcpp::DataFrame &d, size_t col) : modelMatrix(d,col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3];
      hparName = "var." + parName;
   }

   ~modelRanf_cor() {
   }
   


};

#endif /* modelRanf_cor_h */
