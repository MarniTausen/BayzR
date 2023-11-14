//  labeledMatrix.cpp
//

#include "labeledMatrix.h"
#include "rbayzExceptions.h"
#include "nameTools.h"

// add/copy names from Rcpp matrix in the labeledMatrix object.
// Throws errror if rownames not available, auto-fills colnames if colnames not available
// This is now a member function of labeledMatrix, so that the object can call addRowNames
// on itself to get its row or colnames filled.
void labeledMatrix::addRowColNames(Rcpp::NumericMatrix M, std::string & name) {
   rownames = getMatrixNames(M, 1);
   if(rownames.size()==0) {  // rownames empty not allowed
      throw generalRbayzError("No rownames on matrix " + name + "\n");
   }
   colnames = getMatrixNames(M, 2);
   if (rownames.size()==0) { // colnames empty, fill auto colnames
      colnames = generateLabels("col",M.ncol());
   }
}

labeledMatrix::labeledMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(col) {
   // Need to temporarily redo the conversion of the input Robject to
   // Rcpp::NumericMatrix to retrieve row and col names.
   Rcpp::NumericMatrix Rmatrix = Rcpp::as<Rcpp::NumericMatrix>(col);
   this->addRowColNames(Rmatrix, name);
}

void labeledMatrix::initWith(Rcpp::NumericMatrix & M, std::string & name, size_t useCol) {
   simpleMatrix::initWith(M, useCol);
   this->addRowColNames(M, name);
}

