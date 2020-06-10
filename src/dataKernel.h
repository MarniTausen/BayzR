//
//  data_kernel.hpp
//  Data storage class to store a 'kernel' or similarity/relationship matrix.
//  These kernels are now always stored as sets of eigenvectors and are imported
//  via ranf or ran2f functions with arguments V=. The ranf and ran2f already
//  prepare the eigenvector transform.
//  This is an implementation of dataMatrix class filling 'data' with the
//  set of eigen-vectors of the kernel.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataKernel_h
#define dataKernel_h

#include <stdio.h>
#include <Rcpp.h>
#include "dataMatrix.h"

class dataKernel : public dataMatrix {
public:
    dataKernel(Rcpp::RObject col);
    ~dataKernel();
   Rcpp::NumericVector eval;
   size_t nEvalUsed;
   double rrankpct;
};

#endif /* dataKernel_h */
