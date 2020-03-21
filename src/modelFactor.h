//
//  modelTerm_factor.h
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelTerm_factor_h
#define modelTerm_factor_h

#include <Rcpp.h>
#include <cmath>
#include "modelTerm.h"

// Model-term that handles storage of factor data.
// - par vector has size levels of the factor
// - hpar is determined in child classes (not used if fixed, used if random)
// - common working vectors are lhs and rhs vectors
// - define coldata vector, it is IntegerVector for factors
// - common methods are correction, decorrection, and collection of rhs and lhs vectors

class modelTerm_factor : public modelTerm {
   
public:
   
   modelTerm_factor(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      coldata = d[col];
      for (size_t i=0; i<coldata.size(); i++)
         coldata[i] -= 1;
      parLevelNames = coldata.attr("levels");
      par.resize(parLevelNames.size(),0);
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~modelTerm_factor() {
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

   Rcpp::IntegerVector coldata;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations

};

#endif /* modelTerm_factor_h */
