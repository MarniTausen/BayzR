//
//  dataFactor.h
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef dataFactor_h
#define dataFactor_h

#include <Rcpp.h>

class dataFactor {
   
public:
   
   dataFactor(Rcpp::DataFrame &d, size_t col) {
      data = d[col];
      for (size_t i=0; i<data.size(); i++)
         data[i] -= 1;
      labels = coldata.attr("levels");
   }
   
   ~dataFactor() {
   }
   
protected:

   Rcpp::IntegerVector data;
   Rcpp::CharacterVector labels;

};

#endif /* dataFactor_h */
