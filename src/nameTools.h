//
//  BayzR --- nameTools.h
//
//  Tools to retrieve, match names, make indexes, etc.
//  Created by Luc Janss on 02/07/2020.
//

#ifndef nameTools_h
#define nameTools_h

#include <stdio.h>
#include <vector>
#include <string>
#include <Rcpp.h>
#include "dataMatrix.h"
#include "dataFactor.h"

void CharVec2cpp(std::vector<std::string> & labels, Rcpp::CharacterVector templabels);
void getMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat);
void builObsIndex(std::vector<size_t> & obsIndex, dataFactor *F, dataMatrix *M);

#endif /* nameTools_h */
