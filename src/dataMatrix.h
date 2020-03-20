//
//  dataMatrix.h
//  Parent class for storage of (real) matrix, this class has several derived classes,
//  for instance, it can implement kernel matrix, or matrix of covariates.
//  Also nested regressions?
//  The class has generic methods to work on residual corrections and LHS/RHS statistics.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataMatrix_h
#define dataMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "rbayzExceptions.h"

class dataMatrix {
public:
    dataMatrix(Rcpp::RObject col);
    virtual ~dataMatrix();
protected:
   Rcpp::NumericMatrix M;
// double * M;  // want to test difference using low-level C arrays for storage
};

#endif /* dataMatrix_h */
