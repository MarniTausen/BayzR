//
//  modelTerm_2factor.h
//  rbayz
//
//  Created by Luc Janss on 14/01/2020.
//

#ifndef modelTerm_2factor_h
#define modelTerm_2factor_h

#include <Rcpp.h>
#include <cmath>
#include "modelTerm.h"

// Model-term that handles storage of data on 2 factors (used to build classes for interaction fit).
// - par vector size ....
// - hpar is not yet set here, it depends on the actual class if this is used or not
// - common working vectors are lhs and rhs vector
// - define two IntegerVectors for coldata
// - common methods are correction, decorrection, and collection of rhs and lhs vectors

class modelTerm_2factor : public modelTerm {
   
public:
   
   modelTerm_2factor(Rcpp::DataFrame &d, size_t col, Rcpp::RObject &col2) : modelTerm(d, col) {
      col1data = d[col];
      // have to link col2 to col2data vector, but col2 is now RObject
      for (size_t i=0; i<coldata.size(); i++)
         coldata[i] -= 1;
      parLevelNames = coldata.attr("levels");
      par.resize(parLevelNames.size(),0);
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~modelTerm_2factor() {
   }
   
protected:

   void resid_correct() {
      for (size_t obs=0; obs<coldata.size(); obs++)
         resid[obs] -= par[coldata[obs]];
   }
   
   void resid_decorrect() {
      for (size_t obs=0; obs<coldata.size(); obs++)
         resid[obs] += par[coldata[obs]];
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

   Rcpp::IntegerVector col1data, col2data;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* modelTerm_factor_h */
