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
#include <algorithm>

class simpleMatrix {

public:
   
   simpleMatrix() { }
   
   simpleMatrix(size_t nr, size_t nc) {
      doalloc(nr, nc);
      std::fill_n(data0, nr*nc, 0.0l);
   }
   
   simpleMatrix(Rcpp::RObject X) {
      if(Rf_isMatrix(X)) {
         Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(X);
         initWith(M);
      }
      else {
         throw(generalRbayzError("Invalid matrix input in simpleMatrix constructor"));
      }
      
   }
   
   void initWith(Rcpp::NumericMatrix M) {
      // initWith cannot be used on already allocated matrix
      if (nrow> 0 || ncol>0 ) {
         throw(generalRbayzError("Cannot resize matrix in simpleMatrix"));
      }
      doalloc(M.nrow(), M.ncol());
      double *col;  // temporary column pointer
      for(size_t icol=0; icol<M.ncol(); icol++) {
         col=data[icol];
         for(size_t irow=0; irow<M.nrow(); irow++) {
            col[irow] = M(irow,icol);
         }
      }
   }

   ~simpleMatrix() {
      delete[] data;
      delete[] data0;
   }

   double *data0;
   double **data;
   size_t nrow=0,ncol=0;

private:
   // memory alloc may fail, but it will be caught in try-catch in main function
   void doalloc(size_t nr, size_t nc) {
      if (nr <= 0 || nc <= 0) {
         throw(generalRbayzError("Zero or negative sizes in initialisation in simpleMatrix"));
      }
      data0 = new double[nr*nc];    // a block of memory nrows * ncolumns
      data  = new double*[nc];      // a list of column-pointers to index it by column
      for(size_t i=0; i<nc; i++)
         data[i] = data0 + i*nr;
      nrow = nr; ncol = nc;
   }
};

#endif /* simpleMatrix_h */
