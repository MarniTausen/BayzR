//  labeledMatrix.h
//  A simpleMatrix with 'labels' (row and column names).
//  Rownames are required, colnames are optional.
//
//  Created by Luc Janss on 10/05/2021.
//

#ifndef labeledMatrix_h
#define labeledMatrix_h

#include <Rcpp.h>
#include <vector>
#include <string>
#include "rbayzExceptions.h"
#include "simpleMatrix.h"
//#include "nameTools.h"  // strange, including nameTools.h does not work to make
                          // addMatrixNames available here ...

std::vector<std::string> getMatrixNames(Rcpp::NumericMatrix & mat, int dim);

class labeledMatrix : public simpleMatrix {

public:

   labeledMatrix() : simpleMatrix() {  
   }

   // add/copy names from Rcpp matrix in the labeledMatrix object.
   // Throws errror if rownames not available, auto-fills colnames if colnames not available
   addRowColNames(Rcpp::NumericMatrix M, std::string & name) {
      rownames = getMatrixNames(M, 1);
      if(rownames.size()==0) {  // rownames empty not allowed
         throw generalRbayzError("No rownames on matrix " + name + "\n");
      }
      colnames = getMatrixNames(Rmatrix, 2);
      if (rownames.size()==0) { // colnames empty, fill auto colnames
         colnames = generateLabels("col",Rmatrix.ncol());
      }
   }

   labeledMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(col) {
      // Need to temporarily redo the conversion of the input Robject to
      // Rcpp::NumericMatrix to retrieve row and col names.
      Rcpp::NumericMatrix Rmatrix = Rcpp::as<Rcpp::NumericMatrix>(col);
      addRowColNames(Rmatrix, name);
   }

   void initWith(Rcpp::NumericMatrix & M, std::string & name, size_t useCol) {
      simpleMatrix::initWith(M, useCol);
      addRowColNames(M, name);
   }

   ~labeledMatrix() {
   }

   std::vector<std::string> rownames, colnames;
   
};

#endif /* labeledMatrix_h */
