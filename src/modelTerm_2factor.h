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
   
   modelTerm_2factor(Rcpp::DataFrame &d, size_t col) : modelTerm(d, col) {
      Rcpp::RObject thiscol = d[col];
      Rcpp::RObject secondcol = thiscol.attr("factor2");
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[2];
      col1data = d[col];
      col2data = Rcpp::as<Rcpp::NumericVector>(secondcol);
      coldata.push_back(col1data);
      coldata.push_back(col2data);
      for (size_t i=0; i<col1data.size(); i++) col1data[i] -= 1;
      for (size_t i=0; i<col2data.size(); i++) col2data[i] -= 1;
      Rcpp::CharacterVector factor1Names = col1data.attr("levels");
      Rcpp::CharacterVector factor2Names = col2data.attr("levels");
      size_t nLevel1=factor1Names.size(), nLevel2=factor2Names.size();
      // initially set-up interaction coding and solutions for all combinations,
      // factor2 levels nested within factor1 levels. 
      par.resize(nLevel1*nLevel2,0);
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            Rcpp::String s = factor1Names[i];
            s += "%";
            s += factor2Names[j];
            parLevelNames.push_back(s);
         }
      }
      for(size_t i=0; i<col1data.size(); i++) {
         intdata.push_back(col1data[i]*nLevel2 + col2data[i]);
      }
      lhs.resize(parLevelNames.size(),0);
      rhs.resize(parLevelNames.size(),0);
   }
   
   ~modelTerm_2factor() {
   }
   
protected:

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
   
   Rcpp::IntegerVector col1data, col2data, intdata;
   Rcpp::DataFrame coldata;
   std::vector<double> lhs, rhs;

};

#endif /* modelTerm_2factor_h */
