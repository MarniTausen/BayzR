//
//  dataMatrix.h
//  Parent class for storage of (real) matrix, this class has several derived classes,
//  for instance, it can implement kernel matrix, or matrix of covariates.
//  Also nested regressions?
//  The class has generic methods to work on residual corrections and LHS/RHS statistics.
//
//  Data storage class to store a 'kernel' or similarity/relationship matrix.
//  These kernels are now always stored as sets of eigenvectors and are imported
//  via ranf or ran2f functions with arguments V=. The ranf and ran2f already
//  prepare the eigenvector transform.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataMatrix_h
#define dataMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "dataMatrix.h"

class dataMatrix {
public:
    dataMatrix(Rcpp::RObject col);
    virtual ~dataMatrix();
protected:
   Rcpp::NumericMatrix M;
// double * M;  // want to test difference using low-level C arrays for storage
};

#endif /* dataMatrix_h */
