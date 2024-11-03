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
   if (colnames.size()==0) { // colnames empty, fill auto colnames
      colnames = generateLabels("col",M.ncol());
   }
}

// Version restricting the column labels copied to 'useCol' columns, part of constructor
// with 'useCol' limit.
// Proper setting of useCol (>=1 and <= M.ncol) is not tested here, it is assumed that this
// labeling function is only used in labeledMatrix contructor, and then simpleMatrix() constructor
// will throw errors for improper setting of useCol.
void labeledMatrix::addRowColNames(Rcpp::NumericMatrix M, std::string & name, size_t useCol) {
   rownames = getMatrixNames(M, 1);
   if(rownames.size()==0) {  // rownames empty not allowed
      throw generalRbayzError("No rownames on matrix " + name + "\n");
   }
   std::vector<std::string> tempnames = getMatrixNames(M, 2);
   if(tempnames.size()==0) {   // no colnames supplied, auto-generate useCol colnames
      colnames = generateLabels("col",useCol);
   }
   else {                     // copy useCol columns from tempnames in colnames labels
      colnames.resize(useCol);
      for(size_t i=0; i<useCol; i++) colnames[i] = tempnames[i];
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
   this->addRowColNames(M, name, useCol);
}

