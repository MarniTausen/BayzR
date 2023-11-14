//
//  simpleMatrix.h
//  Class to hold a matrix of doubles, quite basic functionality but designed for performance.
//  - all memory allocated as one block
//  - organized column-major so computations go fastest working by column
//  - accessible with double index ->data[col][row] - NOTE col first!
//    or alternatively, ->data[col] is a pointer to the data of column `col`
//  - constructor with sizes specifies number of rows first! (as usual, but opposite to retrieving of elements).
//  - initWith can only be used after using the constructor with no arguments.
// Note: the actual elements need to be accessed as object->data[] (column) or ->data[][] (element),
//    at least the first one can be made more fancy by adding a function for operator[]?
//
//  Created by Luc Janss on 23/10/2020.
//

#ifndef simpleMatrix_h
#define simpleMatrix_h

#include <Rcpp.h>

class simpleMatrix {

public:
   // 'empty' constructor, to be used with one of the initWith() versions   
   simpleMatrix() {   }
   // other constructors
   simpleMatrix(size_t nr, size_t nc);
   simpleMatrix(Rcpp::RObject X);
   
   void initWith(Rcpp::NumericMatrix M, size_t useCol);
   void initWith(Rcpp::NumericMatrix M);

   void swap(simpleMatrix* other);

   ~simpleMatrix();

   double* data0;
   double** data;
   size_t nrow=0,ncol=0;

private:
   void doalloc(size_t nr, size_t nc);

};

#endif /* simpleMatrix_h */
