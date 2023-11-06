//
//  BayzR -- parVector.cpp
//

#include "parVector.h"

// contructor for par-vector with single element where parname is also used for the label
parVector::parVector(std::string parname, double initval) : parVals() {
   parName=parname;
   nelem=1;
   parVals.initWith(1, initval);
   parLabels.push_back(parname);
   val=parVals.data;
   postMean.initWith(nelem,0.0l);
   postSD.initWith(nelem,0.0l);
}

// response model needs a constructor with a vector of values and vector of labels
parVector::parVector(std::string parname, Rcpp::NumericVector initval, Rcpp::CharacterVector& labels)
          : parVals() {
   parName=parname;
   nelem = labels.size();
   parVals.initWith(initval);
   parLabels.resize(labels.size());
   for(Rsize_t i=0; i<labels.size(); i++)
      parLabels[i]=labels[i];
   val=parVals.data;
   postMean.initWith(nelem,0.0l);
   postSD.initWith(nelem,0.0l);
}

// many other model objects can initialize from a single scalar value and labels, the size
// needed is taken from labels size.
parVector::parVector(std::string parname, double initval, Rcpp::CharacterVector& labels)
          : parVals() {
   parName=parname;
   nelem = labels.size();
   parVals.initWith(nelem, initval);
   parLabels.resize(labels.size());
   for(Rsize_t i=0; i<labels.size(); i++)
      parLabels[i]=labels[i];
   val=parVals.data;
   postMean.initWith(nelem,0.0l);
   postSD.initWith(nelem,0.0l);
}

// nearly the same but labels is a vector<string>
parVector::parVector(std::string parname, double initval, std::vector<std::string>& labels)
  : parVals() {
   parName=parname;
   nelem = labels.size();
   parVals.initWith(nelem, initval);
   parLabels.resize(labels.size());
   for(size_t i=0; i<labels.size(); i++)
      parLabels[i]=labels[i];
   val=parVals.data;
   postMean.initWith(nelem,0.0l);
   postSD.initWith(nelem,0.0l);
}

parVector::collecStats() {
  for(size_t i=0; i<nelem; i++) {
    postMean.data[i] += par->val[i];
    postSD[i] += (par->val[i])*(par->val[i]);
  }
  count_collect_stats++;
}

// 
parVector::parInfo() {

}

// Function to write (part of) parVector (for debugging purposes) - seems that Rcout accepts this.
// operator<< overload must be defined as a non-member   
std::ostream& operator<<(std::ostream& os, const parVector& p)
{
    os << p.parName << "[" << p.nelem << "] ";
    size_t loopsize = (p.nelem<5)? p.nelem : 5;
    for(size_t i=0; i<loopsize; i++) os << p.val[i] << " ";
    if(p.nelem >5 ) os << "...";
    os << "\n";
    return os;
}

#endif /* parVector_h */
