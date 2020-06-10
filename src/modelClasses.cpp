//
//  BayzR -- modelClasses.cpp
//
//  Created by Luc Janss on 24/01/2020.
//

#include <Rcpp.h>
#include <cmath>
#include "modelResp.h"
#include "modelMean.h"
#include "modelLiab.h"
#include "modelFactor.h"
#include "modelFixf.h"
#include "modelRanf.h"
#include "modelFreg.h"
#include "modelMatrix.h"











/* old code from ranf_cor for comparison
 void resid_correct(size_t col) {
    for (size_t obs=0; obs < coldata.size(); obs++)
       resid[obs] -= par[col] * matrixdata(coldata(obs),col);
 }
 
 void resid_decorrect(size_t col) {
    for (size_t obs=0; obs < coldata.size(); obs++)
       resid[obs] += par[col] * matrixdata(coldata(obs),col);
 }

 void collect_lhs_rhs(size_t col) {
    lhs = 0.0; rhs=0.0;
    size_t rowlevel;
    for (size_t obs=0; obs < coldata.size(); obs++) {
       rowlevel = coldata(obs);
       rhs += matrixdata(rowlevel,col) * residPrec[obs] * resid[obs];
       lhs += matrixdata(rowlevel,col) * matrixdata(rowlevel,col) * residPrec[obs];
    }
 }

 */
