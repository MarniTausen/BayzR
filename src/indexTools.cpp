//
//  BayzR --- indexTools.cpp
//
//  Created by Luc Janss on 02/07/2020.
//

#include "indexTools.h"
#include "rbayzExceptions.h"

// Comparison function to search for a string in a vector of pair<string, number>
inline bool compString2Pair(const std::pair<std::string, size_t> & p, const std::string & s)
{
    return p.first < s;
}

void builObsIndex(std::vector<size_t> & obsIndex, dataFactor *F, labeledMatrix *M) {
   int errors=0;
//   Rcpp::Rcout << "Going to resize obsIndex to " << F->data.size() << "\n";
   obsIndex.resize(F->data.nelem,0);
   // Build a sorted list of the matrix rownames paired with matrix entry-rows
   std::vector< std::pair<std::string, size_t> > matLabelsSorted;
   for(size_t i=0; i< M->rownames.size(); i++)
      matLabelsSorted.push_back(std::make_pair(M->rownames[i], i));
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
