//
//  BayzR --- indexTools.h
//
//  Tools to match names and make indexes - e.g. between factor and matrix
//  Created by Luc Janss on 02/07/2020.
//
#ifndef indexTools_h
#define indexTools_h

#include <vector>
#include "labeledMatrix.h"
#include "dataFactor.h"

void builObsIndex(std::vector<size_t> & obsIndex, dataFactor *F, labeledMatrix *M);

#endif /* indexTools_h */
