//
//  modelTerm_ran2f_2cor.hpp
//  rbayz
//
//  Created by Luc Janss on 18/01/2020.
//

#ifndef modelTerm_ran2f_2cor_h
#define modelTerm_ran2f_2cor_h

#include <Rcpp.h>
#include "modelTerm_2factor.h"
#include "priorClasses.h"

// Model-term for interaction between two random factors both with covariance-structures.
// The parent class has already set up the interaction coding, here need to add storage
// of relationship matrices and all correct, decorect and lhs-rhs computations need to be modified
// (should be more like ran_cor).
// Note: ran_cor for 1 random variable with correlation derived from realmat class which
// retrieved and stored the matrix. Here we need two matrices and there is no
// class to store 2 matrices ... therefore need to do it here by copying code from realmat ...

class modelTerm_ran2f_2cor : public modelTerm_2factor {

public:

   modelTerm_ran2f_2cor(Rcpp::DataFrame &d, size_t col)  : modelTerm_2factor(d, col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3] + "." + names[4];
      hparName = "var." + parName;
      if (col1data.hasAttribute("evectors")) {
         matrix1data = Rcpp::as<Rcpp::NumericMatrix>(col1data.attr("evectors"));
      }
      else {
         throw(generalRbayzError(std::string("In ran2f: no matrix data attached to variable 1")));
      }
      if (col2data.hasAttribute("evectors")) {
         matrix2data = Rcpp::as<Rcpp::NumericMatrix>(col2data.attr("evectors"));
      }
      else {
         throw(generalRbayzError(std::string("In ran2f: no matrix data attached to variable 2")));
      }
      if(matrix1data.nrow() != matrix1data.ncol()) {
         throw(generalRbayzError(std::string("In ran2f(...,V1=v1) matrix v1 is not square")));
      }
      if(matrix2data.nrow() != matrix2data.ncol()) {
         throw(generalRbayzError(std::string("In ran2f(...,V2=v2) matrix v2 is not square")));
      }
      size_t nLevel1=factor1Names.size(), nLevel2=factor2Names.size();
      if (nLevel1 != matrix1data.nrow()) {
         throw(generalRbayzError(std::string("In ran2f(F1,...,V1=v1) number of levels in F1 do not match size of v1")));
      }
      if (nLevel2 != matrix2data.nrow()) {
         throw(generalRbayzError(std::string("In ran2f(F1,F2,...,V2=v2) number of levels in F2 do not match size of v2")));
      }
      eval1 = Rcpp::as<Rcpp::NumericVector>(col1data.attr("evalues"));
      for (size_t i=0, nPosEval1=0;  i<eval1.size() && eval1[i] >= 0.001; i++) nPosEval1++;
      eval2 = Rcpp::as<Rcpp::NumericVector>(col2data.attr("evalues"));
      // info on the evector columns of the interaction matrix. Every enty in the
      // combination intcol1, intcol2 and evalint tells what column of matrix1data (intcol1)
      // and column of matrix2data (intcol2) it is based on and the eigenvalue (evalint).
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            if (eval1[i]*eval2[j] > 0.001) {
               evalint.push_back(eval1[i]*eval2[j]);
               intcol1.push_back(i);
               intcol2.push_back(j);
            }
         }
      }
      Rcpp::Rcout << "ran2f with two V-matrices models interaction VxV with rank " << evalint.size() << std::endl;
      if(evalint.size()==0)
         throw generalRbayzError("Zero rank in ran2f for VxV; all eigenvalues are below tolerance?");
      workcol.resize(nLevel1*nLevel2,0);
   }

   ~modelTerm_ran2f_2cor() {
   }

   void sample() {
      size_t nLevel1=factor1Names.size(), nLevel2=factor2Names.size();
      size_t matrix1col, matrix2col;
      for(size_t k=0; k<evalint.size(); k++) {
         matrix1col = intcol1[k];  // column k of interaction matrix is combination of these
         matrix2col = intcol2[k];  // two columns of the two input relationship matrices
         for(size_t i=0; i<nLevel1; i++) {
            for(size_t j=0; j<nLevel2; j++) {
               workcol[i*nLevel2+j] = matrix1data(i,matrix1col) * matrix2data[j,matrix2col];
            }
         }
         // the workcol still needs to be multiplied over the factor, from this
         // point all should work the same as in realmat....
         
         resid_decorrect();
         collect_lhs_rhs();
         for(size_t k=0; k<par.size(); k++) {
            // random effect: add 1/hpar[0] in lhs
            par[k] = R::rnorm( (rhs[k]/(lhs[k]+(1/hpar[0]))), sqrt(1.0/(lhs[k]+(1/hpar[0]))));
         }
         resid_correct();
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<par.size(); k++)
         ssq += par[k]*par[k];
      hpar[0] = gprior.samplevar(ssq,par.size());
   }

private:

   // this is code from ran2f that collects all lhs and rhs at once, probably
   // needs to be changed to collect one level at-a-time like ran_cor.
   void collect_lhs_rhs() {
      size_t k;
      for(k=0; k<par.size(); k++) {
         rhs[k] = 0.0;
         lhs[k] = 0.0;
      }
      for (size_t obs=0; obs<intdata.size(); obs++) {
         k=intdata[obs];
         rhs[k] += residPrec[obs] * resid[obs];
         lhs[k] += residPrec[obs];
      }
   }
   
   Rcpp::NumericMatrix matrix1data, matrix2data;
   Rcpp::NumericVector eval1, eval2;
   std::vector<double> evalint, workcol;
   std::vector<size_t> intcol1,intcol2;
   size_t nPosEval1, nPosEval2, nPosEvali;


};

#endif /* modelTerm_ran2f_2cor_h */
