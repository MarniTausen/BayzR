//
//  dataMatrix.h
//  Parent class for storage of (real) matrix. It is quite empty but defines
//  that all derived classes have a NumericMatrix 'data' that stores a matrix.
//  Class modelMatrix has matching computational operation on this 'data'.
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
   dataMatrix(Rcpp::RObject col) { }
   virtual ~dataMatrix();
   Rcpp::NumericMatrix data;
   // double * data;  // want to test difference using low-level C arrays for storage
};

#endif /* dataMatrix_h */
