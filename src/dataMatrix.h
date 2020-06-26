//
//  dataMatrix.h
//  Class for storage of matrix input data: a kernel/similarity matrix, or matrix of covariates.
//  The different types of input are loaded based on attributes found in the column data.
//  Class modelMatrix has matching computational operations on this 'data'.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataMatrix_h
#define dataMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "rbayzExceptions.h"

class dataMatrix {

public:
   dataMatrix(Rcpp::RObject col) {
      if (col.hasAttribute("evectors")) {
         double rrankpct;
         data = Rcpp::as<Rcpp::NumericMatrix>(col.attr("evectors"));
         if (col.hasAttribute("evalues")) {
            weights = Rcpp::as<Rcpp::NumericVector>(col.attr("evalues"));
            // the weigts should actually be the inverse eigen-values, but it is dangerous
            // to start making the inverses already here because some eigenvalues can be zero ...
         }
         else
            throw(generalRbayzError("Eigen-decomp storage corrupted - no evalues"));
         if(data.nrow() != data.ncol())
            throw(generalRbayzError("Eigen-decomp has non-square e-vectors matrix"));
         if(weights.size() != data.ncol())
            throw(generalRbayzError("Eigen-decomp vectors and values size do not match"));
         if(col.hasAttribute("rrankpct"))
            rrankpct = col.attr("rrankpct");
         else
            rrankpct = 99;
         double sumeval = 0.0l;
         for (size_t i = 0; i<weights.size() && weights[i] > 0; i++) sumeval += weights[i];
         double eval_cutoff = rrankpct * sumeval / 100.0l;
         nColUsed = 0;
         sumeval = 0.0l;
         while (sumeval < eval_cutoff) sumeval += weights[nColUsed++];
         // this message can be moved to the modelTerm routine; and it must be different for ran2f_2cor.
         Rcpp::Rcout << "In ranf with V rrankpct=" << rrankpct << " uses " << nColUsed << "eigenvectors\n";
      }
      else {
         Rcpp::Rcout << "Matrix input not (yet) supported on column xxx\n";
      }
   }

   ~dataMatrix() {
   }

   Rcpp::NumericMatrix data;
   Rcpp::NumericVector weights;
   size_t nColUsed;

   // double * data;  // want to test difference using low-level C arrays for storage
};

#endif /* dataMatrix_h */
