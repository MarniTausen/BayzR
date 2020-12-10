//
//  BayzR --- nameTools.cpp
//
//  Created by Luc Janss on 02/07/2020.
//

#include "nameTools.h"
#include "rbayzExceptions.h"
#include <algorithm>
#include <map>

//#include "dataMatrix.h"
//#include "dataFactor.h"

// Transform R CharacterVector to the c++ equivalent vector<string>
// Note: the C++ vector should be empty, the algorithm uses push_back.

void CharVec2cpp(std::vector<std::string> & CppStrings, Rcpp::CharacterVector RStrings) {
   std::string s;
   for(size_t i=0; i< RStrings.size(); i++) {
      s = Rcpp::as<std::string>(RStrings[i]);
      CppStrings.push_back(s);
   }
}

// Comparison function to search for a string in a vector of pair<string, number>
inline bool compString2Pair(const std::pair<std::string, size_t> & p, const std::string & s)
{
    return p.first < s;
}

// getMatrixNames: retrieve row or col-names (use dim=1 or 2 as in apply()) from an Rcpp
//     NumericMatrix, returns NULL vector when dimnames not present or the requested dim is empty.
// addMatrixNames: fills matrix row or col-names (dim=1 or 2) in an c++ vector<string>,
//     return value 0 for success, 1 for errors (missing names).
// Note: addMatrixNames uses push_back, typical use is to pass an empty vector<string> as 'names'
// argument and it will be filled.

Rcpp::CharacterVector getMatrixNames(Rcpp::NumericMatrix & mat, int dim) {
   if (mat.hasAttribute("dimnames")) {
      Rcpp::List dimnames = Rcpp::as<Rcpp::List>(mat.attr("dimnames"));
      if(dim==1) return dimnames[0];
      else return dimnames[1];
   }
   else return R_NilValue;
}

int addMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat, int dim) {
   Rcpp::CharacterVector matNames = getMatrixNames(mat, dim);
   if(matNames.isNULL()) return 1;
   else {
      CharVec2cpp(names, matNames);
      return 0;
   }
}

// Find 'name' in the column-names of a data frame.
// The column names are not sorted, so applying a simple sequential search, but it will
// get slow if the data frame is large (e.g. with large covariate data in it!).
// The c++ find will do the same and does not return the index, so I just made a simple
// search myself that goes over the vector of names.
// Could be improved by making a sorted version of the column names and using binary_search,
// but ideally then also storing and re-using the sorted version for all look-ups.
// Then the sorted names must be prepared in the main functinon....?
int findDataColumn(Rcpp::DataFrame d, std::string name) {
   std::vector<std::string> colnames;
   CharVec2cpp(colnames, d.names());
   int col;
   for(col=0; col<colnames.size(); col++) {
      if(colnames[col]==name) break;
   }
   if(col==colnames.size())
      return -1;
   else
      return col;
}

void builObsIndex(std::vector<size_t> & obsIndex, dataFactor *F, dataMatrix *M) {
   int errors=0;
//   Rcpp::Rcout << "Going to resize obsIndex to " << F->data.size() << "\n";
   obsIndex.resize(F->data.nelem,0);
   // Build a sorted list of the matrix labels paired with matrix entry-rows
   std::vector< std::pair<std::string, size_t> > matLabelsSorted;
   for(size_t i=0; i< M->labels.size(); i++)
      matLabelsSorted.push_back(std::make_pair(M->labels[i], i));
   std::sort(matLabelsSorted.begin(),matLabelsSorted.end()); // shoud be enough
   // Make the index that links observations to matrix rows
   std::string s;
   std::vector< std::pair<std::string, size_t> >::iterator it;
   for(size_t i=0; i<F->data.nelem; i++) {
      s = F->labels[F->data[i]];
      if ( (it = std::lower_bound(matLabelsSorted.begin(), matLabelsSorted.end(), s, compString2Pair)) != matLabelsSorted.end()) {
         obsIndex[i]=it->second;
      }
      else {
         // error
         Rcpp::Rcout << "ID " << s << "is not available in similarity/kernel matrix\n";
         errors++;
      }
   }
   if(errors)
      throw (generalRbayzError("Some IDs not found in similarity/kernel matrix"));
}
