//  RcppFactor
//

#include "RcppFactor.h"

RcppFactor (Rcpp::RObject col) {

   if(Rf_isFactor(col))  : levels() {
      Rcpp::IntegerVector tempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(tempvec);
      levels.initWith(tempvec);
      Rcpp::CharacterVector templabels = col.attr("levels");
      CharVec2cpp(labels, templabels);
      for(size_t row=0; row < unsigned(tempvec.size()); row++) {
         levels[row] -= 1;
      }
      // If there was "NA" in the data, add it as extra (last) level
      if (Rcpp::sum(missing) > 0) {
         labels.push_back("NA");
         size_t last_level = labels.size()-1;
         for(size_t row=0; row < unsigned(tempvec.size()); row++) {
            if(missing[row]) data[row] = last_level;
         }
      }
   }
   else if (Rcpp::is<Rcpp::NumericVector>(col) && !Rf_isMatrix(col)) {

   }
   else if (Rcpp::is<Rcpp::CharacterVector>(col) && !Rf_isMatrix(col)) {

   }
   else if(Rcpp::is<Rcpp::LogicalVector>(col)  && !Rf_isMatrix(col)) {

   }
   else {
      // error: not convertable to factor
   }

}

