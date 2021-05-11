//  dataMatrix.h
//  Class for storage of matrix of covariates (e.g. input in rr() model).
//  Derives from labeledMatrix so it has rownames (and optionally colnames), and
//  makes column centering of the input data.
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

class dataMatrix : public labeledMatrix {

public:
   dataMatrix(Rcpp::RObject col, std::string & name) : labeledMatrix(col, name) {
      // column-center the matrix data
      double * datacol;
      size_t i,j;
      double sum;
      for(i=0; i<ncol; i++) {
         datacol = data[i];
         sum=0.0l;
         for(j=0; j<nrow; j++) sum+=datacol[j];
         sum /= double(nrow);
         for(j=0; j<nrow; j++) datacol[j] -= sum;
      }
   }

   ~dataMatrix() {
   }

};

#endif /* dataMatrix_h */
