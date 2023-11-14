//
//  simpleMatrix.cpp
//

#include "simpleMatrix.h"
#include <algorithm>
#include "rbayzExceptions.h"

   
simpleMatrix::simpleMatrix(size_t nr, size_t nc) {
   doalloc(nr, nc);
   std::fill_n(data0, nr*nc, 0.0l);
}

simpleMatrix::simpleMatrix(Rcpp::RObject X) {
   if(Rf_isMatrix(X)) {
      Rcpp::NumericMatrix M = Rcpp::as<Rcpp::NumericMatrix>(X);
      initWith(M);
   }
   else {
      throw(generalRbayzError("Invalid matrix input in simpleMatrix constructor"));
   }
   
}

// initialise (allocate and fill) matrix from existing Rcpp matrix,
// useCol sets the number of columns to use/copy.
// There is also an initWith without useCol to copy the whole Rcpp matrix.
void simpleMatrix::initWith(Rcpp::NumericMatrix M, size_t useCol) {
   if (nrow> 0 || ncol>0 ) {    // oops this matrix is already allocated
      throw(generalRbayzError("Cannot resize matrix in simpleMatrix"));
   }
   if (useCol > unsigned(M.ncol())) {
      throw(generalRbayzError("useCol is larger than actual ncol in simpleMatrix"));
   }
   doalloc(M.nrow(), useCol);
   double *col;  // temporary column pointer
   for(size_t icol=0; icol<useCol; icol++) {
      col=data[icol];
      for(size_t irow=0; irow< unsigned(M.nrow()); irow++) {
         col[irow] = M(irow,icol);
      }
   }
}
void simpleMatrix::initWith(Rcpp::NumericMatrix M) {
   initWith(M, M.ncol());
}
// Swap contents of two matrices: the contents of this-> (object itself) are
// swapped with content of matrix pointed to by other->. 
void simpleMatrix::swap(simpleMatrix* other) {
   double** olddata = this->data;
   double* olddata0 = this->data0;
   size_t oldnrow   = this->nrow;
   size_t oldncol   = this->ncol;
   this->data  = other->data;
   this->data0 = other->data0;
   this->nrow  = other->nrow;
   this->ncol  = other->ncol;
   other->data = olddata;
   other->data0 = olddata0;
   other->nrow  = oldnrow;
   other->ncol  = oldncol;
}
~simpleMatrix::simpleMatrix() {
   delete[] data;
   delete[] data0;
}
// memory alloc. This may fail, but it will be caught in try-catch in main function
void simpleMatrix::doalloc(size_t nr, size_t nc) {
   if (nr <= 0 || nc <= 0) {
      throw(generalRbayzError("Zero or negative sizes in initialisation in simpleMatrix"));
   }
   data0 = new double[nr*nc];    // a block of memory nrows * ncolumns
   data  = new double*[nc];      // a list of column-pointers to index it by column
   for(size_t i=0; i<nc; i++)
      data[i] = data0 + i*nr;
   nrow = nr; ncol = nc;
}
