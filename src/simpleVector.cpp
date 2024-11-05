//
//  simpleVector.cpp
//

#include "simpleVector.h"
#include <algorithm>
#include "rbayzExceptions.h"

simpleIntVector::simpleIntVector(size_t n) {
   doalloc(n);
   std::fill_n(data, n, 0);
}
simpleIntVector::~simpleIntVector() {
   if (nelem>0) {
      delete[] data;
   }
}

int& simpleIntVector::operator[](size_t i) {
   return(data[i]);
}
void simpleIntVector::initWith(Rcpp::IntegerVector v) {
   doalloc(v.size());
   for(size_t i=0; i<nelem; i++)
      data[i] = v[i];
}
void simpleIntVector::initWith(size_t n, int initvalue) {
   doalloc(n);
   for(size_t i=0; i<nelem; i++)
      data[i] = initvalue;
}
void simpleIntVector::doalloc(size_t n) {
   if (n <= 0) {
      throw(generalRbayzError("Zero or negative size in initialisation in simpleVector"));
   }
   if (nelem >0) {  // already allocated!
      throw(generalRbayzError("Attempting realloc in simpleVector not supported"));
   }
   data  = new int[n];
   nelem = n;
}

// simpleDblVector functions
// -------------------------
simpleDblVector::simpleDblVector(size_t n) {
   doalloc(n);
   std::fill_n(data, n, 0.0l);
}
simpleDblVector::~simpleDblVector() {
   if (nelem>0) {
      delete[] data;
   }
}

double& simpleDblVector::operator[](size_t i) {
   return(data[i]);
}

void simpleDblVector::initWith(Rcpp::NumericVector v, size_t useElem) {
   if (useElem > unsigned(v.size())) {
      throw(generalRbayzError("useElem is larger than actual nelem in simpleDblVector"));
   }
   doalloc(useElem);
   for(size_t i=0; i<useElem; i++)
      data[i] = v[i];
}

void simpleDblVector::initWith(Rcpp::NumericVector v) {
   initWith(v, v.size());
}

void simpleDblVector::initWith(size_t n, double initvalue) {
   doalloc(n);
   for(size_t i=0; i<nelem; i++) data[i] = initvalue;
}

void simpleDblVector::initWith(simpleDblVector & X) {
   doalloc(X.nelem);
   for(size_t i=0; i<nelem; i++) data[i] = X.data[i];
}

// [ToDo] make a wrap feature for simpleVector? - it can avoid copying.
// A safe way to guarantee that the vector being wrapped does not disappear (making
// the wrapper fail), is to put a "lock" in the one being wrapped that avoids clean-up?
// But can it work when cleanup is needed to get the order right that first the wrapper
// is removed and the one being wrapped is released? 

void simpleDblVector::swap(simpleDblVector* other) {
   double* olddata = this->data;
   size_t oldnelem   = this->nelem;
   this->data  = other->data;
   this->nelem  = other->nelem;
   other->data = olddata;
   other->nelem  = oldnelem;
}

void simpleDblVector::doalloc(size_t n) {
   if (n <= 0) {
      throw(generalRbayzError("Zero or negative size in initialisation in simpleVector"));
   }
   if (nelem >0) {  // already allocated!
      throw(generalRbayzError("Attempting realloc in simpleVector not supported"));
   }
   data  = new double[n];
   nelem = n;
}

/*
I tried to merge the code of simpleIntVector and simpleDblVector by using templates,
and I got all code units compiling, but in the end the linker still complained and it
was not clear how to resolve that.

The class with templates looked like:
template <typename T> class simpleVector {
public:
   simpleVector() : nelem(0) { }
   simpleVector(size_t n);
   ~simpleVector();
   T& operator[](size_t i);    // retrieve and set value using []
   void initWith(size_t n, T initvalue);
   template <typename T2> void initWith(T2 v);
   template <typename T2> void initWithPart(T2 v, size_t useElem);
   void swap(simpleVector* other);
   T *data;
   size_t nelem;
private:
   void doalloc(size_t n);
};

Some of the methods looked like:
template <typename T> T& simpleVector<T>::operator[](size_t i) {
   return(data[i]);
}
template <typename T> template<typename T2> void simpleVector<T>::initWith(T2 v) {
   doalloc(v.size());
   for(size_t i=0; i<nelem; i++)
      data[i] = v[i];
}

Some of the linker errors:
undefined reference to `void simpleVector<int>::initWith<Rcpp::Vector<13, Rcpp::PreserveStorage> >(Rcpp::Vector<13, Rcpp::PreserveStorage>)'
undefined reference to `simpleVector<int>::operator[](unsigned int)
undefined reference to `simpleVector<int>::~simpleVector()
*/
