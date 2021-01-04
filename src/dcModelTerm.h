//
//  dcModelTerm.h
//  A class to describe a DeConstructed Model-Term - where all elements of a model-term
//  are taken apart and stored in separate slots.
//
//  Created by Luc Janss on 29/12/2020.

#ifndef dcModelTerm_h
#define dcModelTerm_h

#include <Rcpp.h>
#include <vector>
#include <string>

class dcModelTerm {
public:
   dcModelTerm(std::string modelTerm, Rcpp::DataFrame &d);
   ~dcModelTerm() { }
   std::vector<std::string> variableNames;
   std::vector<Rcpp::RObject> variableObjects;
   std::vector<int> variableTypes;
   enum varianceType {notgiven, matrices, linmod};
   enum hierarchType {nohierarch, index, genuine};
   std::vector<std::string> varianceNames;
   std::string varianceModel;
   std::string priorModel;
   std::string hierarchModel;
};

#endif /* dcModelTerm */
