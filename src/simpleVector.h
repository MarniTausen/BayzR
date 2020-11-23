//
//  simpleVector.h
//  Wrapper around a basic C-style vector (pointer to block of memory), but with
//  C++ allocation using new. Now with version for Integer data (simpleIntVector) and
//  preparing also one for doubles (SimpleDblVector).
//  Developing now two uses:
//   - if at declaration size is known, can use constructor with size and it will initialise with 0
//   - if vector is class member and size is not known at construction, constructor with no arguments
//     can be used, it will not allocated the vector yet. Then use initWith() to both allocated and
//     initialise it with inpput from an Rcpp vector.
//  To access the vector elements use .data[i] or ->data[i], size is stored as member variable nelem.
//
//  Created by Luc Janss on 23/10/2020.
//

#ifndef simpleIntVector_h
#define simpleIntVector_h

#include <Rcpp.h>
#include <algorithm>

class simpleIntVector {

public:
   
   simpleIntVector() : nelem(0) {
   }
   
   simpleIntVector(size_t n) {
      if (nelem <= 0) {
         throw(generalRbayzError("Zero or negative size in initialisation in simpleVector"));
      }
      nelem = n;
      data  = new int[n];
      std::fill_n(data, n, 0);
   }

   ~simpleIntVector() {
      if (nelem>0) {
         delete[] data;
      }
   }
   
   int& operator[](size_t i) {
      return(data[i]);
   }
   
   void initWith(Rcpp::IntegerVector v) {
      if(nelem>0) {
         throw(generalRbayzError("Cannot initialise already allocated vector in simpleVector"));
         // can be extended to re-allocate
      }
      size_t n = v.size();
      data  = new int[n];
      nelem = n;
      for(size_t i=0; i<nelem; i++)
         data[i] = v[i];
   }

   int *data;
   size_t nelem;
   
};



#endif /* simpleIntVector_h */
