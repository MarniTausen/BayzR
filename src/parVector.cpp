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
   postVar.initWith(nelem,0.0l);
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
   postVar.initWith(nelem,0.0l);
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
   postVar.initWith(nelem,0.0l);
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
   postVar.initWith(nelem,0.0l);
}

// Update cumulative means and variances
void parVector::collecStats() {
   double olddev, newdev;
   count_collect_stats++;
   double n = (double)count_collect_stats;
   if (count_collect_stats==1) {                  // only update mean
      for(size_t i=0; i<nelem; i++) {
         olddev = par->val[i] - postMean.data[i];
         postMean.data[i] += olddev/n;
      }
   }
   else {                                         // can update mean and var
      for(size_t i=0; i<nelem; i++) {
         olddev = par->val[i] - postMean.data[i]; // deviation current value with old mean
         postMean.data[i] += olddev/n;
         newdev = par->val[i] - postMean.data[i]; // deviation with updated mean
         postVar[i] += (olddev*newdev)/(n-1.0l);
      }
   }
}

// 
void parVector::setParInfo(std::string fname, std::string vars, std::string varstruct) {
   modelFunction=fname;
   variableList=vars;
   varianceStruct=varstruct;
}

// Function to write (part of) parVector (for debugging purposes) - Rcout accepts this fine.
// operator<< overload must be defined as a non-member ...
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
