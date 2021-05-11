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

int addMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat, int dim);

class labeledMatrix : public simpleMatrix {

public:

   labeledMatrix() : simpleMatrix() {  }

   labeledMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(col) {
      // I need to temporarily redo the conversion of the input Robject to
      // Rcpp::NumericMatrix to retrieve row and col names.
      Rcpp::NumericMatrix tempdata = Rcpp::as<Rcpp::NumericMatrix>(col);
      if (addMatrixNames(rownames, tempdata, 1) >0) {
         throw generalRbayzError("No rownames on matrix " + name + "\n");
      }
      addMatrixNames(colnames, tempdata, 2); // no throw here, colnames are optional
   }

   void initWith(Rcpp::NumericMatrix M, size_t useCol) {
      simpleMatrix::initWith(M, useCol);
      if (addMatrixNames(rownames, M, 1) >0) {
         throw generalRbayzError("No rownames on matrix xxx\n");  // no name available here!
      }
      addMatrixNames(colnames, M, 2); // no throw here, colnames are optional

   }

   ~labeledMatrix() {
   }

   std::vector<std::string> rownames, colnames;
   
};

#endif /* labeledMatrix_h */
