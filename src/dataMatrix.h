//
//  dataMatrix.h
//  Class for storage of matrix input data (matrix of covariates).
//  This is a slightly specialised (derived) version of simpleMatrix: it also has
//  names for the rows, and weights.
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
#include <string>
#include <vector>
#include "simpleMatrix.h"
#include "simpleVector.h"
//#include "nameTools.h"  // strange, including nameTools.h does not work to make
                          // addMatrixNames available here ...

int addMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat, int dim);

class dataMatrix : public simpleMatrix {

public:
   dataMatrix(Rcpp::RObject col) : simpleMatrix(col) {
      // constuctor of simpleMatrix has already done the matrix data, here add
      // weights and labels, and column-center the covariates.
      Rcpp::NumericMatrix tempdata = Rcpp::as<Rcpp::NumericMatrix>(col);
      V.initWith(tempdata.ncol(), 1.0l);
      weights = V.data;  // the weights pointer in a dataMatrix object is the 'data' in V
      if (addMatrixNames(rownames, tempdata, 1) >0) {
         throw generalRbayzError("No rownames on matrix ... \n");
      }
      addMatrixNames(colnames, tempdata, 2); // colnames are allowed to be missing
      double * datacol;
      size_t i,j;
      double sum;
      for(i=0; i<ncol; i++) {
         datacol = data[i];
         sum=0.0l;
         for(j=0; j<nrow; j++) sum+=datacol[j];
         for(j=0; j<nrow; j++) datacol[j] -= sum;
      }
      /*  this part needs to be reviewed, it needs to move to variance-input
      if (col.hasAttribute("evectors")) {
         throw(generalRbayzError("Evector input as matrix data temporarily disables\n"));
         double rrankpct;
         data = Rcpp::as<Rcpp::NumericMatrix>(col.attr("evectors"));
         std::vector<std::string> att = data.attributeNames();
         Rcpp::Rcout << "Handling evector matrix with attributes: ";
         for(size_t i=0; i<att.size(); i++) Rcpp::Rcout << " " << att[i];
         Rcpp::Rcout << "\n";
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
         Rcpp::Rcout << "In ranf with V rrankpct=" << rrankpct << " uses " << nColUsed << " eigenvectors\n";
         getMatrixNames(labels, data);  // can be moved outside if-else when other matrix forms are also coded
      }
      else  {
       */
   }

   ~dataMatrix() {
   }

   simpleDblVector V;
   double *weights;
   std::vector<std::string> rownames, colnames;
   
};

#endif /* dataMatrix_h */
