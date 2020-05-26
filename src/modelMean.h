//
//  BayzR --- modelMean.h
//
// Computational class for fitting mean (intercept)
// - this one has no data!
// - par is size 1 to hold mean
// - hpar is not used
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef modelMean_h
#define modelMean_h

#include <Rcpp.h>


class modelMean : public modelBase {
   
public:
   
   modelMean(Rcpp::DataFrame &d, size_t col) : modelBase(d, col) {
      par.resize(1,0);
      parName = "mean";
   }

   ~modelMean() {
   }

   void sample();

private:
};

#endif /* modelMean_h */
