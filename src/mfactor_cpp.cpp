//
//  mfactor_cpp.cpp
//
//  Created by Luc Janss on 29/01/2020.
//
//  Function to create a 'multifactor' from multiple columns input data;
//  a multi-factor is an integer matrix with factor coding across columns: there is one
//  'levels' attribute connected for all the unique levels across columns.
//  Input columns can be: factor, integer, character.

#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>

// [[Rcpp::export]]
Rcpp::IntegerVector mfactor_cpp(Rcpp::List fact_list) {

   // Check lengths and class of input columns. Input columns should be factor, integer or character,
   // but it is enough to check integer or character (factor is also integer).
   size_t nr_cols = fact_list.length();
   size_t nr_rows=0, l=0;
   for (size_t col=0; col<nr_cols; col++) {
      if (Rcpp::is<Rcpp::IntegerVector>(fact_list[col])) {
         Rcpp::IntegerVector v = fact_list[col];
         l = v.length();
      }
      else if (Rcpp::is<Rcpp::CharacterVector>(fact_list[col])) {
         Rcpp::CharacterVector v = fact_list[col];
         l = v.length();
      }
      else {  // not right class for this column
         std::string s = "Column ";
         s += std::to_string(col+1);
         s += " is not factor, integer or character; try as.character()?\n";
         Rcpp::stop(s);
      }
      if(col==0) nr_rows=l;
      else {
         if(l != nr_rows) // unequal lengths
            Rcpp::stop("Input columns do not all have the same length\n");
      }
   }
   
   // Make a C++ vector of strings to collect all column data as c++ text strings
   std::vector<std::string> all_input;
   all_input.reserve(nr_rows*nr_cols);
   for (size_t col=0; col<nr_cols; col++) {
      if(Rf_isFactor(fact_list[col])) {
         Rcpp::IntegerVector v = fact_list[col];
         Rcpp::CharacterVector levelNames = v.attr("levels");
         Rcpp::String s;
         for(size_t row=0; row<nr_rows; row++) {
            s = levelNames[v[row]-1];
            all_input.push_back(s);
         }
      }
      else if (Rcpp::is<Rcpp::IntegerVector>(fact_list[col])) {
         Rcpp::IntegerVector v = fact_list[col];
         std::string s;
         for(size_t row=0; row<nr_rows; row++) {
            s = std::to_string(v[row]);
            all_input.push_back(s);
         }
      }
      else if (Rcpp::is<Rcpp::CharacterVector>(fact_list[col])) {
         Rcpp::CharacterVector v = fact_list[col];
         Rcpp::String s;
         for(size_t row=0; row<nr_rows; row++) {
            s = v[row];
            all_input.push_back(s);
         }
      }
   }
   
   // Now build a map of the unique levels across all input columns
   std::map<std::string, int> merged_levels;
   std::map<std::string, int>::iterator p;
   for(size_t i=0; i<all_input.size(); i++) {
      p = merged_levels.lower_bound(all_input[i]);
      if (p==merged_levels.end() || all_input[i] < p->first) // conditions for not in map
         merged_levels.insert(p,std::make_pair(all_input[i],0));
   }
   
   // Code the merged levels in the map, starting from 1 to make it compatible with R
   l=1;
   for(p=merged_levels.begin(); p != merged_levels.end(); p++) {
      p->second = l++;
   }
   
   // Then create an R integer vector with coded levels; it is transformed later
   // into a matrix to make the final result.
   Rcpp::IntegerVector mfactor(nr_cols*nr_rows);
   for(size_t i=0; i<all_input.size(); i++) {
      p = merged_levels.find(all_input[i]);
      mfactor[i] = p->second;
   }

   // Copy the levels stored in the map to R CharacterVector
   Rcpp::CharacterVector new_levels(merged_levels.size());
   l=0;
   for(p=merged_levels.begin(); p != merged_levels.end(); p++) {
      new_levels[l++] = p->first;
   }

   mfactor.attr("dim") = Rcpp::Dimension(nr_rows,nr_cols);
   mfactor.attr("levels") = new_levels;
   mfactor.attr("class") = "mfactor";
   return mfactor;
   
   
}



