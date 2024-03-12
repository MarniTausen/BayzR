//  simpleFactor
//

#include "simpleFactor.h"
#include "nameTools.h"
#include "rbayzExceptions.h"

simpleFactor::simpleFactor (Rcpp::RObject col, std::string inp_name)  : simpleIntVector() {

/* simpleFactor can take different types of input and convert it to a 'factor', but
   processing and handling NAs differs somewhat dependent on input type. I always
   keep NA as the last level.
   - if input is R factor, it is directly convertable to integer vector that can be
     copied in the 'data' using unitWith, but NAs become a large negative number and are
     repaired later.
   - if input is an R character vector, converting to strings would convert NAs to string "NA",
     but I avoid that by building a map with the input strings that are not NA. Then later
     NA is added as last level, and it's level-codes are inserted as the last 'lev' counter value.
   - if input is an R integer vector I build a map of <int,int> so the sorting is nicer (as
     strings it would give the ugly sorting with "10" before "2"); treatment of NAs is like for
     a character vector.
   - logical vector input can be directly converted to interger and copied in 'data' using
     initWith, NAs will become large negative numbers and needs NA treatment like the first case.
     The IntegerVector will have 0 for false, 1 for true, so it is immediate in C base-0 coding.
*/
   name = inp_name;
   if(Rf_isFactor(col))  {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      initWith(Rtempvec);  // init of 'data' vector
      Rcpp::CharacterVector templabels = col.attr("levels");
      CharVec2cpp(labels, templabels);
      for(size_t row=0; row < unsigned(Rtempvec.size()); row++) {
         data[row] -= 1;
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
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      initWith(Rtempvec);
      std::map<int, int> unique_levels;
      for(int i=0; i<Rtempvec.size(); i++) {  // build map with non-NA values
         if(!missing[i]) unique_levels[Rtempvec[i]];
      }
      std::map<int, int>::iterator p;
      lev=0;        // Code the merged levels in the map
      for(p=unique_levels.begin(); p != unique_levels.end(); p++) p->second = lev++;
      initWith(Rtempvec.size(), 0.0l);
      for(int i=0; i<Rtempvec.size(); i++) {  // code the data
         if(missing[i])
            data[i] = lev;
         else {
            p = unique_levels.find(Rtempvec[i]);
            data[i] = p->second;
         }
      }
      // fill labels vector
      labels.reserve(lev);
      for(p=unique_levels.begin(); p != unique_levels.end(); p++)
         labels.push_back(std::to_string(p->first));
      if(sum(missing)>0) labels_push_back("NA");
   }
   else if (Rcpp::is<Rcpp::CharacterVector>(col) && !Rf_isMatrix(col)) {
      Rcpp::CharacterVector Rtempvec = Rcpp::as<Rcpp::CharacterVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      std::vector<std::string> Ctempvec;
      CharVec2cpp(Ctempvec, Rtempvec); 
      std::map<std::string, int> unique_levels;
      for(size_t i=0; i<Ctempvec.size(); i++) {  // build map with non-NA values
         if(!missing[i]) unique_levels[Ctempvec[i]];
      }
      std::map<std::string, int>::iterator p;
      lev=0;        // Code the merged levels in the map
      for(p=unique_levels.begin(); p != unique_levels.end(); p++) p->second = lev++;
      initWith(Ctempvec.size(), 0.0l);
      for(size_t i=0; i<Ctempvec.size(); i++) {  // code the data
         if(missing[i])
            data[i] = lev;
         else {
            p = unique_levels.find(Ctempvec[i]);
            data[i] = p->second;
         }
      }
      // fill labels vector
      labels.reserve(lev);
      for(p=unique_levels.begin(); p != unique_levels.end(); p++) labels.push_back(p->first);
      if(sum(missing)>0) labels_push_back("NA");
   }
   else if(Rcpp::is<Rcpp::LogicalVector>(col)  && !Rf_isMatrix(col)) {
      Rcpp::IntegerVector Rtempvec = Rcpp::as<Rcpp::IntegerVector>(col);
      Rcpp::LogicalVector missing = Rcpp::is_na(Rtempvec);
      initWith(Rtempvec);  // init of 'data' vector
      labels.push_back("FALSE");
      labels.push_back("TRUE");
      if (Rcpp::sum(missing) > 0) {
         labels.push_back("NA");
         for(size_t row=0; row < unsigned(tempvec.size()); row++) {
            if(missing[row]) data[row] = 2;
         }
      }
   }
   else {
      throw generalRbayzError("Variable/data column is not convertable to a factor: " + name);
   }
}

// Convert stored factor data back to the 'full' vector of strings.
std::vector<std::string> simpleFactor::back2vecstring() {
   std:vector<std::string> result(nelem);
   for(size_t i=0; i<nelem; i++) {
      result[i] = labels[data[i]];
   }
   return result;
}