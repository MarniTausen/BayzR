//
//  modelTerm_realmat.h
//  rbayz
//
//  Created by Luc Janss on 30/08/2019.
//

#ifndef modelTerm_realmat_h
#define modelTerm_realmat_h

#include <Rcpp.h>
#include <cmath>
#include "dataMatrix.h"
#include "modelBase.h"

// Model-term where 'coldata' is a matrix of real numbers (Rcpp NumericMatrix, C++ double), it is
// a parent class of model-terms for random regression and also the ran_cor (that works like
// a random regression model on eigenvectors). Here:
// - par vector is size number of columns of the input matrix (will hold regression coefficients)
// - hpar is not set here, it depends on the actual class if this is used or not ... <- check
// - common working vectors are lhs and rhs vector
// - the matrix that should be attached to coldata should be passed as argument (reference) in the
//   constructor
// - common methods are correction, decorrection, and collection of rhs and lhs vectors

class modelMatrix : public modelBase {
   
public:
   
   // The constructor can work in different 'modes' for receiving and storing matrix-data.
   // flag=1: (to be added) for covariate data
   // flag=2: for data from ran_cor where column data is a factor and there is a correlation
   //         matrix that comes in as eigen-vectors
   // .... ?
   
   modelTerm_realmat(Rcpp::DataFrame &d, size_t col, int flag) : modelTerm(d, col) {
      if (flag==2) {
         coldata = d[col];
         for (size_t i=0; i<coldata.size(); i++)
            coldata[i] -= 1;
         parLevelNames = coldata.attr("levels");
         par.resize(parLevelNames.size(),0);
         Rcpp::RObject thiscol = d[col];
         matrixdata = Rcpp::as<Rcpp::NumericMatrix>(thiscol.attr("evectors"));
         if(matrixdata.nrow() != matrixdata.ncol()) {
            throw(generalRbayzError(std::string("In ranf(...,V=v) matrix v is not square")));
         }
          // this check can more to ran_cor
         if (parLevelNames.size() != matrixdata.nrow()) {
            throw(generalRbayzError(std::string("In ranf(f,V=v) number of levels in f do not match size of v")));
         }
         
         /* This is code for change to the 'U-tilde' model
         Rcpp::NumericVector d = Rcpp::as<Rcpp::NumericVector>(thiscol.attr("evalues"));
         double evalsqrt;
         for(size_t col=0; col<coldata.ncol(); col++) {
            if (d[col] >= 0.001) {
               evalsqrt = sqrt(d[col]);
            }
            else {
               evalsqrt=0.0;
            }
            for(size_t row=0; row<coldata.nrow(); row++)
               coldata(row,col) = coldata(row,col)*evalsqrt;
         }
         */
      }
   }
   
   ~modelTerm_realmat() {
   }
   
protected:

   void resid_correct(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] -= par[col] * matrixdata(coldata(obs),col);
   }
   
   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] += par[col] * matrixdata(coldata(obs),col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      size_t rowlevel;
      for (size_t obs=0; obs < coldata.size(); obs++) {
         rowlevel = coldata(obs);
         rhs += matrixdata(rowlevel,col) * residPrec[obs] * resid[obs];
         lhs += matrixdata(rowlevel,col) * matrixdata(rowlevel,col) * residPrec[obs];
      }
   }

   dataMatrix M;
   Rcpp::IntegerVector coldata;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   
};

#endif /* modelTerm_realmat_h */
