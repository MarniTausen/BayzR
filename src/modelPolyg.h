//
//  BayzR --- modelTerm.hpp
//
// Model-term for polygenic effect - not yet ready, needs a pedigree data class?
//
// The rbayz wrapper function makes this column a character vector, and stores pedigree
// as an extra attribute in the column data.
// - the coldata is CharacterVector
// - par vector is size of pedigree (but should be unique)
// - hpar is size 1 and used for genetic variance

//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelPolyg_h
#define modelPolyg_h

#include <Rcpp.h>

class modelPolyg : public modelBase {

public:

   modelPolyg(parsedModelTerm & modeldescr, modelBase * rmod)
         : modelBase(modeldescr, rmod) {
      coldata = d[col];
      // need to process the column data, simple character list? Needs index or searching?
      // + need to extract and process pedigree
      par = new parVector(modeldescr, 0.0l); // not right yet
      // variance needs to be added
   }

   ~modelPolyg() {
   }

   void sample() {
      
   }

private:
   Rcpp::IntegerVector coldata;  // probably need to change, I want to prepare column for polyg() as character

};

#endif /* modelPolyg_h */
