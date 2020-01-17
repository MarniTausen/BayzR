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

// Model-term where 'coldata' is a matrix of real numbers (Rcpp NumericMatrix, C++ double), it is
// a parent class of model-terms for random regression and also the ran_cor (that works like
// a random regression model on eigenvectors). Here:
// - par vector is size number of columns of the input matrix (will hold regression coefficients)
// - hpar is not set here, it depends on the actual class if this is used or not ... <- check
// - common working vectors are lhs and rhs vector
// - the matrix that should be attached to coldata should be passed as argument (reference) in the
//   constructor
// - common methods are correction, decorrection, and collection of rhs and lhs vectors

class modelTerm_realmat : public modelTerm {
   
public:
   
   modelTerm_realmat(Rcpp::DataFrame &d, size_t col, int flag) : modelTerm(d, col) {
      Rcpp::RObject thiscol = d[col];
      if (flag==2) {  // this indicates the matrix data is stored as evectors and evalues
         coldata = Rcpp::as<Rcpp::NumericMatrix>(thiscol.attr("evectors"));
         /* This is code if would like to change to the 'U-tilde' model
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
      parLevelNames = colnames(coldata);
      par.resize(parLevelNames.size(),0);
   }
   
   ~modelTerm_realmat() {
   }
   
protected:

   void resid_correct(size_t col) {
      for (size_t obs=0; obs < coldata.nrow(); obs++)
         resid[obs] -= par[col] * coldata(obs,col);
   }
   
   void resid_decorrect(size_t col) {
      for (size_t obs=0; obs < coldata.nrow(); obs++)
         resid[obs] += par[col] * coldata(obs,col);
   }

   void collect_lhs_rhs(size_t col) {
      lhs = 0.0; rhs=0.0;
      for (size_t obs=0; obs < coldata.nrow(); obs++) {
         rhs += residPrec[obs] * resid[obs] * coldata(obs,col);
         lhs += residPrec[obs] * residPrec[obs];
      }
   }

   Rcpp::NumericMatrix coldata;
   double lhs, rhs;          // lhs, rhs will be scalar here (per iteration)
   std::vector<double> fit;
   
};

#endif /* modelTerm_realmat_h */
