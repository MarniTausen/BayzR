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
#include <math.h>
#include "rbayzExceptions.h"
#include "simpleMatrix.h"

class dataMatrix : public labeledMatrix {

public:
   dataMatrix(Rcpp::RObject col, std::string & name) : labeledMatrix(col, name) {
      // column-center the matrix data and fill missings with column mean.
      // (and after centering the column means are zero for all columns).
      double * datacol;
      size_t i,j, nobs;
      double sum;
      for(i=0; i<ncol; i++) {
         datacol = data[i];
         sum=0.0l; nobs=0;
         for(j=0; j<nrow; j++) {
            if (! isnan(datacol[j]) ) { // cannot catch it as NA_REAL after is has been copied
               sum+=datacol[j];         // in the simpleMatrix structure ... maybe can be improved?
               nobs++;
            }
         }
         sum /= double(nobs);
         for(j=0; j<nrow; j++) {
            if (isnan(datacol[j]))
               datacol[j] = 0.0l;
            else 
               datacol[j] -= sum;
         }
      }
   }

   ~dataMatrix() {
   }

};

#endif /* dataMatrix_h */
