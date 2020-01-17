//
//  modelTerm_realvec.h
//  rbayz
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelTerm_realvec_h
#define modelTerm_realvec_h

#include <Rcpp.h>
#include <cmath>
#include "modelTerm.h"

// Model-term where 'coldata' is a vector of real numbers (Rcpp NumericVector, C++ double), it is
// a parent class of the model-term for fixed regression. Here:
// - par vector is size 1 (but can also choose to set the size in fixreg class)
// - hpar is not used
// - common working vectors are lhs and rhs vector
// - common methods are correction, decorrection, and collection of rhs and lhs vectors

class modelTerm_realvec : public modelTerm {
   
public:
   
   modelTerm_realvec(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      coldata = Rcpp::as<Rcpp::NumericVector>(d[col]);
//      parLevelNames = colnames(coldata);
      par.resize(1,0);
   }
   
   ~modelTerm_realvec() {
   }
   
protected:

   void resid_correct() {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] -= par[0] * coldata[obs];
   }
   
   void resid_decorrect() {
      for (size_t obs=0; obs < coldata.size(); obs++)
         resid[obs] += par[0] * coldata[obs];
   }

   void collect_lhs_rhs() {
      lhs = 0.0; rhs=0.0;
      for (size_t obs=0; obs < coldata.size(); obs++) {
         rhs += residPrec[obs] * resid[obs] * coldata[obs];
         lhs += coldata[obs] * residPrec[obs] * coldata[obs];
      }
   }

   Rcpp::NumericVector coldata;
   double lhs, rhs;
   std::vector<double> fit;  // not used yet
   
};

#endif /* modelTerm_realvec_h */
