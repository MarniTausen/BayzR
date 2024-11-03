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
#include "simpleMatrix.h"

class labeledMatrix : public simpleMatrix {

public:

   // 'empty' constructor, to be used with initWith()
   labeledMatrix() : simpleMatrix() {  
   }

   // constructor with an R matrix
   labeledMatrix(Rcpp::RObject col, std::string & name);
   void addRowColNames(Rcpp::NumericMatrix M, std::string & name);
   void addRowColNames(Rcpp::NumericMatrix M, std::string & name, size_t useCol);

   // this is hiding the initWith methods from simpleMatrix class, and a labeledMatrix
   // can only be initialised with this one?
   void initWith(Rcpp::NumericMatrix & M, std::string & name, size_t useCol);

   ~labeledMatrix() {
   }

   std::vector<std::string> rownames, colnames;
   
};

#endif /* labeledMatrix_h */
