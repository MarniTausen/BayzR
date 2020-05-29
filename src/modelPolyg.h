//
//  modelTerm.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef modelTerm_polyg_h
#define modelTerm_polyg_h

#include <Rcpp.h>

// Model-term for polygenic effect.
// The rbayz wrapper function makes this column a character vector, and stores pedigree
// as an extra attribute in the column data.
// - the coldata is CharacterVector
// - par vector is size of pedigree (but should be unique)
// - hpar is size 1 and used for genetic variance

class modelTerm_polyg : public modelTerm {

public:

   modelTerm_polyg(Rcpp::DataFrame &d, size_t column) : modelTerm(d, col) {
      coldata = d[col];
      // need to process the column data, simple character list? Needs index or searching?
      // + need to extract and process pedigree
      hpar.resize(1,0);
      hparName = "var." + parName;
   }

   ~modelTerm_polyg() {
   }

   void sample() {
      
   }

private:
   Rcpp::IntegerVector coldata;  // probably need to change, I want to prepare column for polyg() as character

};

#endif /* modelTerm_polyg_h */
