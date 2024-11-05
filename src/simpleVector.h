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

class simpleIntVector {
public:
   simpleIntVector() { }
   simpleIntVector(size_t n);
   virtual ~simpleIntVector();
   int& operator[](size_t i);    // retrieving and setting value using []
   void initWith(Rcpp::IntegerVector v);
   void initWith(size_t n, int initvalue);
   int *data;
   size_t nelem=0;
private:
   void doalloc(size_t n);
};

class simpleDblVector {
public:
   simpleDblVector() { }
   simpleDblVector(size_t n);
   virtual ~simpleDblVector();
   double& operator[](size_t i);
   void initWith(Rcpp::NumericVector v, size_t useElem);
   void initWith(Rcpp::NumericVector v);
   void initWith(size_t n, double initvalue);
   void initWith(simpleDblVector & X);
   void swap(simpleDblVector* other);
   double *data;
   size_t nelem=0;
private:
   void doalloc(size_t n);
};


#endif /* simpleVector_h */
