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
~simpleIntVector::simpleIntVector() {
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
      throw(generalRbayzError("Cannot resize vector in simpleVector"));
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
~simpleDblVector::simpleDblVector() {
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
   for(size_t i=0; i<nelem; i++)
      data[i] = initvalue;
}
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
      throw(generalRbayzError("Cannot resize vector in simpleVector"));
   }
   data  = new double[n];
   nelem = n;
}
