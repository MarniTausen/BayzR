//
//  modelTerm_ran2f_2cor.hpp
//  rbayz
//
//  Created by Luc Janss on 18/01/2020.
//

#ifndef modelTerm_ran2f_2cor_h
#define modelTerm_ran2f_2cor_h

#include <Rcpp.h>
#include "model2Factor.h"
#include "priorClasses.h"
#include <vector>

// Model-term for interaction between two random factors both with covariance-structures.
// The parent class has already set up the interaction coding, here need to add storage
// of relationship matrices and all correct, decorect and lhs-rhs computations need to be modified
// (should be more like ran_cor).
// Note: ran_cor for 1 random variable with correlation derived from realmat class which
// retrieved and stored the matrix. Here we need two matrices and there is no
// class to store 2 matrices ... therefore need to do it here by copying code from realmat ...

class modelRan2f_2cor : public model2factor {

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
      eval2 = Rcpp::as<Rcpp::NumericVector>(col2data.attr("evalues"));
	  // Note: when computing evals for Kronecker product from eval1*eval2 entries, these are not sorted!
	  // First compute all product eigenvalues, then sort and determine cut-off, then re-fill in unsorted order with only the needed ones.
	  evalint.resize(nLevel1*nLevel2);
      intcol1.resize(nLevel1*nLevel2);
      intcol2.resize(nLevel1*nLevel2);
	  double sumeval = 0.0l;
	  for (size_t i = 0; i<nLevel1; i++) {
		  for (size_t j = 0; j<nLevel2; j++) {
			  if (eval1[i] <= 0.0l || eval2[j] <= 0.0l)
				  evalint[i*nLevel2 + j] = 0.0l;
			  else
			      evalint[i*nLevel2 + j] = eval1[i] * eval2[j];
			  sumeval += evalint[i*nLevel2 + j];
		  }
	  }
	  std::sort(evalint.begin(), evalint.end(), std::greater<double>());
	  rrankpct = col1data.attr("rrankpct");
	  double eval_sum_cutoff = rrankpct * sumeval / 100.0l;  // this is cut-off on cumulative/sum of eval
	  sumeval = 0.0l;
	  unsigned long nEvalUsed=0;
	  while ((sumeval += evalint[nEvalUsed]) < eval_sum_cutoff) nEvalUsed++;
	  double eval_min_cutoff = evalint[nEvalUsed];  // this is cut-off evalue below which not to use
	  nEvalUsed = 0;
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            if ( !(eval1[i] <= 0.0l || eval2[j] <= 0.0l) && (eval1[i]*eval2[j] >= eval_min_cutoff) ) {
               evalint[nEvalUsed]=eval1[i]*eval2[j];
               intcol1[nEvalUsed]=i;
               intcol2[nEvalUsed]=j;
               nEvalUsed++;
            }
         }
      }
      // Instead of screen this could go to some kind of 'notes' buffer
      Rcpp::Rcout << "ran2f with two V-matrices rrankpct=" << rrankpct << " uses " << nEvalUsed << " eigenvectors\n";
      if(nEvalUsed==0)
         throw generalRbayzError("Zero rank in ran2f for VxV; all eigenvalues are below tolerance?");
      workcol.resize(nLevel1*nLevel2,0);
   }

   ~modelTerm_ran2f_2cor() {
   }

   void sample() {
      size_t nLevel1=factor1Names.size(), nLevel2=factor2Names.size();
      size_t matrix1col, matrix2col, rowlevel;
      double lhsl, rhsl; // local scalar version, there is also vector lhs and rhs in the object
      for(size_t k=0; k<evalint.size(); k++) {
         matrix1col = intcol1[k];  // column k of interaction matrix is combination of these
         matrix2col = intcol2[k];  // two columns of the two input relationship matrices
         for(size_t i=0; i<nLevel1; i++) {
            for(size_t j=0; j<nLevel2; j++) {
               workcol[i*nLevel2+j] = matrix1data(i,matrix1col) * matrix2data(j,matrix2col);
            }
         }
//         Rcpp::Rcout << "Workcol " << k << ": ";
//         for(size_t i=0; i<workcol.size(); i++)
//            Rcpp::Rcout << workcol[i] << " ";
//         Rcpp::Rcout << std::endl;
         // the workcol still needs to be multiplied over the factor, from this
         // point all should work the same as in realmat....
         for (size_t obs=0; obs < intdata.size(); obs++) {
            resid[obs] += par[k] * workcol[intdata[obs]];
         }
         lhsl = 0.0; rhsl=0.0;
         for (size_t obs=0; obs < intdata.size(); obs++) {
            rowlevel = intdata[obs];
            rhsl += workcol[rowlevel] * residPrec[obs] * resid[obs];
            lhsl += workcol[rowlevel] * workcol[rowlevel] * residPrec[obs];
         }
         lhsl = lhsl + (1.0/(evalint[k]*hpar[0]));  // lhs with variance added
         par[k] = R::rnorm( (rhsl/lhsl), sqrt(1.0/lhsl));
         for (size_t obs=0; obs < intdata.size(); obs++)
            resid[obs] -= par[k] * workcol[intdata[obs]];
      }
      // update hyper-par (variance) using SSQ of random effects
      double ssq=0.0;
      for(size_t k=0; k<evalint.size(); k++)
         ssq += par[k]*par[k]/evalint[k];
      hpar[0] = gprior.samplevar(ssq,evalint.size());
   }

private:
   
   Rcpp::NumericMatrix matrix1data, matrix2data;
   Rcpp::NumericVector eval1, eval2;
   std::vector<double> evalint, workcol;
   std::vector<size_t> intcol1,intcol2;
   double rrankpct;


};

#endif /* modelTerm_ran2f_2cor_h */
