//
//  dataMatrix.h
//  Class for storage of matrix of covariates (e.g. input in rr() model).
//  Derives from simpleMatrix and using simpleMatrix(Robject) constructor so that the
//  input is already stored in parent class constructor.
//  Here add row and optionally column names and center the input data.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataMatrix_h
#define dataMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include <string>
#include "rbayzExceptions.h"
#include "simpleMatrix.h"
//#include "nameTools.h"  // strange, including nameTools.h does not work to make
                          // addMatrixNames available here ...

int addMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat, int dim);

class dataMatrix : public simpleMatrix {

public:
   dataMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(col) {
      // I need to temporarily redo the conversion of the input Robject to
      // Rcpp::NumericMatrix to retrieve row and col names.
      Rcpp::NumericMatrix tempdata = Rcpp::as<Rcpp::NumericMatrix>(col);
      if (addMatrixNames(rownames, tempdata, 1) >0) {
         throw generalRbayzError("No rownames on matrix " + name + "\n");
      }
      addMatrixNames(colnames, tempdata, 2); // no throw here, colnames are optional
      double * datacol;
      size_t i,j;
      double sum;
      for(i=0; i<ncol; i++) {
         datacol = data[i];
         sum=0.0l;
         for(j=0; j<nrow; j++) sum+=datacol[j];
         sum =/ double(nrow);
         for(j=0; j<nrow; j++) datacol[j] -= sum;
      }
   }

   ~dataMatrix() {
   }

   std::vector<std::string> rownames, colnames;
   
};

#endif /* dataMatrix_h */
