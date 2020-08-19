//
// rbayz ---- modelRan2f_2cor.hpp
// Model-class for interaction between two random factors both with covariance-structures.
// Note: does not derive from modelMatrix (like Ranf_cor), but from 2Factor, the computations
// and set-up are quite different from the modelMatrix methods because here we have 2 matrices.
// Note2: with correlation matrices, the levels of the factor variable are no longer the levels
// in the data, but become the names (and order) of the names on the matrices, there can be more
// levels in the matrices than in the data.
//
// Created by Luc Janss on 18/01/2020.
//

#ifndef modelRan2f_2cor_h
#define modelRan2f_2cor_h

#include <Rcpp.h>
#include "model2Factor.h"
#include "dataMatrix.h"
#include "nameTools.h"
#include <vector>

class modelRan2f_2cor : public model2Factor {

public:

   modelRan2f_2cor(Rcpp::DataFrame &d, size_t col)  : model2Factor(d, col) {
      hpar.resize(1,1);
      std::vector<std::string> names = parseColNames(d,col);
      parName = parName + "." + names[3] + "." + names[4];
      hparName = "var." + parName;
      Rcpp::RObject col1 = d[col];
      Rcpp::RObject col2 = col1.attr("factor2");
      M1 = new dataMatrix(col1);
      M2 = new dataMatrix(col2);
      size_t nLevel1=M1->data.nrow();           // note levels are now the matrix size
      size_t nLevel2=M2->data.nrow();
      par.resize(nLevel1*nLevel2,0);            // redo these sizes, model2Factor constructor
      parLabels.resize(nLevel1*nLevel2,"");     // may have set them wrong
      for(size_t i=0; i<nLevel1; i++) {         // ... and redo these labels based on matrix
         for(size_t j=0; j<nLevel2; j++) {      // names and order
            Rcpp::String s = M1->labels[i];
            s += "%";
            s += M2->labels[j];
            parLabels[i*nLevel2+j]=s;
         }
      }
      // compute all eigenvalues (and their sum) of the interaction (Kronecker product) matrix
 	   evalint.resize(nLevel1*nLevel2);
      double sumeval = 0.0l;
	   for (size_t i = 0; i<nLevel1; i++) {
		  for (size_t j = 0; j<nLevel2; j++) {
			  if (M1->weights[i] <= 0.0l || M2->weights[j] <= 0.0l)
				  evalint[i*nLevel2 + j] = 0.0l;
			  else
			      evalint[i*nLevel2 + j] = M1->weights[i] * M2->weights[j];
			  sumeval += evalint[i*nLevel2 + j];
		  }
      }
      // I want to determine a cut-off below which eval not to use, but the list is not sorted,
      // so I need to sort first, then run a cumulative sum and determine at what eval I have the
      // desired portion of the total sum.
	   std::sort(evalint.begin(), evalint.end(), std::greater<double>());
	   rrankpct = col1.attr("rrankpct");
	   double eval_sum_cutoff = rrankpct * sumeval / 100.0l;  // this is cut-off on cumulative/sum of eval
	   sumeval = 0.0l;
	   unsigned long nEvalUsed=0;
	   while ((sumeval += evalint[nEvalUsed]) < eval_sum_cutoff) nEvalUsed++;
	   double eval_min_cutoff = evalint[nEvalUsed];  // this is cut-off evalue below which not to use
      // now can build information on the interaction-evectors to use: from which M1 and M2 column
      // each derives (intcol1 and intcol2 gives entries), and the eigenvalue for that combination.
      // Combinations falling below the cut-off are skipped and do not get in the list.
      evalint.resize(nEvalUsed+1);
      intcol1.resize(nEvalUsed+1);
      intcol2.resize(nEvalUsed+1);
      nEvalUsed = 0;
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            if ( !(M1->weights[i] <= 0.0l || M2->weights[j] <= 0.0l)
                    && (M1->weights[i]*M2->weights[j] >= eval_min_cutoff) ) {
               evalint[nEvalUsed]=M1->weights[i] * M2->weights[j];
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
      // also redo the 'intdata' vector from 2Factor class
      // Build observation index from observation indices to the M1 and M2 matrix from temporary
      // indexes that link observations to each matrix
      std::vector<size_t> obsIndex1, obsIndex2;  // these are local and will go out of scope
      builObsIndex(obsIndex1,F1,M1);             // when constructor finishes
      builObsIndex(obsIndex2,F2,M2);
      for(size_t obs=0; obs<F1->data.size(); obs++) {
         intdata[obs] = obsIndex1[obs]*nLevel2 + obsIndex2[obs];
      }
      workcol.resize(nLevel1*nLevel2,0);
   }

   ~modelRan2f_2cor() {
   }

   // A function to make or set one interaction-evector column in the 'workcol' member variable.
   // By putting this in a function I can modify the strategy how to store or compute
   // this interaction-evector, currently it is done 'on the fly' from the two matrices.
   void makeIntCol(size_t col) {
      size_t matrix1col = intcol1[col];  // column 'col' of interaction matrix is combination of these
      size_t matrix2col = intcol2[col];  // two columns of the two input relationship matrices
      size_t nLevel1=M1->labels.size();
      size_t nLevel2=M2->labels.size();
      for(size_t i=0; i<nLevel1; i++) {
         for(size_t j=0; j<nLevel2; j++) {
            workcol[i*nLevel2+j] = M1->data(i,matrix1col) * M2->data(j,matrix2col);
         }
      }
   }
   
   void sample() {
      size_t rowlevel;
      double lhsl, rhsl; // local scalar version, there is also vector lhs and rhs in the object
      for(size_t k=0; k<evalint.size(); k++) {
         makeIntCol(k);
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
   
   dataMatrix *M1, * M2;
   std::vector<double> evalint, workcol;
   std::vector<size_t> intcol1,intcol2;
   double rrankpct;

};

#endif /* modelRan2f_2cor_h */
