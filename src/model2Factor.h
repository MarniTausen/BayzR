//
//  rbayz -- model2Factor.h
//  Computational methods to work on two factors.
//
//  Created by Luc Janss on 29/05/2020.
//

#ifndef model2Factor_h
#define model2Factor_h

#include <Rcpp.h>
#include <cmath>
#include "modelBase.h"
#include "dataFactor.h"

class model2Factor : public modelBase {
   
public:
   
   model2Factor(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      Rcpp::RObject col1 = d[col];
      Rcpp::RObject col2 = col1.attr("factor2");
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[2];
      F1 = new dataFactor(col1);
      F2 = new dataFactor(col2);
      size_t nLevel1=F1->labels.size();
      size_t nLevel2=F2->labels.size();
      // at the moment this still sets-up interaction coding and solutions for
      // all combinations ...
      par.resize(nLevel1*nLevel2,0);
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            Rcpp::String s = F1->labels[i];
            s += "%";
            s += F2->labels[j];
            parLevelNames.push_back(s);
         }
      }
      for(size_t i=0; i<F1->data.size(); i++) {
         intdata.push_back(F1->data[i]*nLevel2 + F2->data[i]);
      }
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~model2Factor() {
      delete F1;
      delete F2;
   }
   
   // The correction and de-correction methods are general for all child classes.
   // The collect_lhs_rhs() methods are different between the child classes.
   void resid_correct() {
      for (size_t obs=0; obs<intdata.size(); obs++)
         resid[obs] -= par[intdata[obs]];
   }
   
   void resid_decorrect() {
      for (size_t obs=0; obs<intdata.size(); obs++)
         resid[obs] += par[intdata[obs]];
   }

   void collect_lhs_rhs();

   dataFactor *F1, *F2;
   std::vector<double> lhs, rhs;          // working vectors to collect LHS an RHS of equations
   Rcpp::IntegerVector intdata;
};

#endif /* model2Factor_h */
