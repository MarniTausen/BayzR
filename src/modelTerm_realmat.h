//
//  modelTerm_realmat.h
//  rbayz
//
//  Created by Luc Janss on 30/08/2019.
//  Copyright © 2019 Luc Janss. All rights reserved.
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
   
   modelTerm_realmat(Rcpp::DataFrame &d, size_t col, Rcpp::NumericMatrix &m) : modelTerm(d, col) {
      coldata = m;
      parLevelNames = colnames(m);
      par.resize(parLevelNames.size(),0);
//      lhs.resize(parLevelNames.size(),0);
//      rhs.resize(parLevelNames.size(),0);
 */
   }
   
   ~modelTerm_realmat() {
   }
   
protected:

   void resid_correct() {
      for (size_t obs=0; obs<coldata.size(); obs++)
         resid[obs] -= par[coldata[obs]];
   }
   
   void resid_decorrect() {
      for (size_t obs=0; obs<coldata.size(); obs++)
         resid[obs] -= par[coldata[obs]];
   }

   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.size(); k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs<coldata.size(); obs++) {
         k=coldata[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }

   Rcpp::NumericMatrix coldata;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations
// this one can probably also have a fit vector
   
};

#endif /* modelTerm_realmat_h */
