//
//  simpleVector.h
//  Wrapper around a basic C-style vector (pointer to block of memory), but with
//  C++ allocation using new. Version for Integer data (simpleIntVector, alternative for Rcpp
//  IntegerVector) and a one for doubles (SimpleDblVector, alternative for Rcpp NumericVector).
//  Developing now two uses:
//   - if at declaration size is known, can use constructor with size and it will initialise with 0
//   - if vector is class member and size is not known at construction, constructor with no arguments
//     can be used, it will not allocated the vector yet. Then use initWith() to both allocated and
//     initialise it with inpput from an Rcpp vector.
//  To access the vector elements use .data[i] or ->data[i], size is stored as member variable nelem.
//
//  Created by Luc Janss on 23/10/2020.
//

#ifndef simpleVector_h
#define simpleVector_h

#include <Rcpp.h>
#include <algorithm>

class simpleIntVector {

public:
   
   simpleIntVector() : nelem(0) { }
   
   simpleIntVector(size_t n) {
      doalloc(n);
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

/*   void operator=(int v) {
      for(size_t i=0; i<nelem; i++)
         data[i] = v;
   }
*/

   void initWith(Rcpp::IntegerVector v) {
      doalloc(v.size());
      for(size_t i=0; i<nelem; i++)
         data[i] = v[i];
   }

   void initWith(size_t n, int initvalue) {
      doalloc(n);
      for(size_t i=0; i<nelem; i++)
         data[i] = initvalue;
   }

   int *data;
   size_t nelem;
   
private:
   void doalloc(size_t n) {
      if (n <= 0) {
         throw(generalRbayzError("Zero or negative size in initialisation in simpleVector"));
      }
      if (nelem >0) {  // already allocated!
         throw(generalRbayzError("Cannot resize vector in simpleVector"));
      }
      data  = new int[n];
      nelem = n;
   }
   
};

class simpleDblVector {

public:
   
   simpleDblVector() : nelem(0) { }
   
   simpleDblVector(size_t n) {
      doalloc(n);
      std::fill_n(data, n, 0.0l);
   }

   ~simpleDblVector() {
      if (nelem>0) {
         delete[] data;
      }
   }
   
   double& operator[](size_t i) {
      return(data[i]);
   }

/*   void operator=(double v) {
      for(size_t i=0; i<nelem; i++)
         data[i] = v;
   }
*/
   
   void initWith(Rcpp::NumericVector v) {
      doalloc(v.size());
      for(size_t i=0; i<nelem; i++)
         data[i] = v[i];
   }

   void initWith(size_t n, double initvalue) {
      doalloc(n);
      for(size_t i=0; i<nelem; i++)
         data[i] = initvalue;
   }

   double *data;
   size_t nelem;

private:
   void doalloc(size_t n) {
      if (n <= 0) {
         throw(generalRbayzError("Zero or negative size in initialisation in simpleVector"));
      }
      if (nelem >0) {  // already allocated!
         throw(generalRbayzError("Cannot resize vector in simpleVector"));
      }
      data  = new double[n];
      nelem = n;
   }

};


#endif /* simpleVector_h */
