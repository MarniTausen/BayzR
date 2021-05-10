//  kernelMatrix.h
//  Storage of kernel as eigen-decomposition.
//  Derives from simpleMatrix using the simpleMatrix() contructor, so that
//  no matrix data is stored yet from the parent constructor.
//
//  Created by Luc Janss on 06/05/2021.
//

#ifndef kernelMatrix_h
#define kernelMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "rbayzExceptions.h"
#include <string>
#include <vector>
#include "simpleMatrix.h"
#include "simpleVector.h"

class kernelMatrix : public simpleMatrix {

public:

   kernelMatrix(Rcpp::RObject col, std::string & name) : labeledMatrix(), weights() {
      kerneldata = Rcpp::as<Rcpp::NumericMatrix>(col);
   	Rcpp::Function eig("eigen");
	   Rcpp::List eigdecomp;
   	try {
	   	eigdecomp = eig(kerneldata);
   	}
   	catch(std::exception &err) {
         throw(generalRbayzError("An error occurred running eigen(): "+err.what()));
   	}
   	Rcpp::NumericVector eigvalues = eigdecomp["values"];
	   Rcpp::NumericMatrix eigvectors = eigdecomp["vectors"];
      // Determine how many evectors to keep
      double rrankpct=90;          // rrankpct not coming correctly from the modelterm now
      size_t nColUsed = 0;
      double sumeval = 0.0l;                   // only summing the positive ones!
      for (size_t i = 0; i<eigvalues.size() && eigvalues[i] > 0; i++) sumeval += eigvalues[i];
      double eval_cutoff = rrankpct * sumeval / 100.0l;
      sumeval = 0.0l;
      while (sumeval < eval_cutoff) sumeval += eigvalues[nColUsed++];
      // this message can be moved to the modelTerm routine; and it must be different for ran2f_2cor.
      Rcpp::Rcout << "For kernel " << name << " with rrankpct=" << rrankpct << " using "
                  << nColUsed << " eigenvectors\n";
      this->initWith(eigvectors, nColUsed);
      weights.initWith(eigvalues, nColUsed);
   }

   ~kernelMatrix() {
   }

   // Add a kernel (make the kronecker product) to the stored kernel in the object.
   void addKernel(kernelMatrix* K2) {
      // make list (and sort it) of all interaction evalues
 	   std::vector<double> evalint(this->ncol*K2->ncol, 0.0l);
      double sumeval = 0.0l;
	   for (size_t i = 0; i<this->ncol; i++) {
		  for (size_t j = 0; j<K2->ncol; j++) {
           evalint[i*K2->ncol + j] = this->weights[i] * K2->weights[j];
			  sumeval += evalint[i*K2->ncol + j];
		  }
      }
	   std::sort(evalint.begin(), evalint.end(), std::greater<double>());
      // Determine cut-off eval to reach rrankpct cumulative sum
      double rrankpct=90;          // rrankpct not coming correctly from the modelterm now
	   double eval_sum_cutoff = rrankpct * sumeval / 100.0l;  // this is cut-off on cumulative eval
	   sumeval = 0.0l;
	   unsigned long nEvalUsed=0;
	   while ((sumeval += evalint[nEvalUsed]) < eval_sum_cutoff) nEvalUsed++;
	   double eval_min_cutoff = evalint[nEvalUsed];  // this is cut-off evalue to get rrankpct cumulative
      nEvalUsed++;
      // Set-up a temporary  matrix and vector to make the evectors and evalues for the interaction kernel.
      // For the moment I keep setting up for all combinations. 
      size_t nLevel1 = this->nrow;
      size_t nLevel2 = K2->nrow;
      simpleMatrix tempEvecs(nLevel1*nLevel2, nEvalUsed);
      simpleDblVector tempEvals(nEvalUsed);
      size_t k = 0;
      for(size_t i=0; i<this->ncol; i++) {
         for(size_t j=0; j<K2->ncol; j++) {
            if ( this->weights[i]*K2->weights[j] >= eval_min_cutoff ) {
               tempEvals[k] = this->weights[i] * K2->weights[j];
               for(size_t rowi=0; rowi<nLevel1; rowi++) {
                  for(size_t rowj=0; rowj<nLevel2; rowj++) {
                     tempEvecs[i*nLevel2+j,k] = this->data(rowi,i) * M2->data(rowj,j);
                  }
               }
               k++;
            }
         }
      }
      // Make labels for the new combined matrix
      std::vector<std::string> tempLabels;
      tempLabels.reserve(nLevel1*nLevel2);
      for(size_t rowi=0; rowi<nLevel1; rowi++) {
         for(size_t rowj=0; rowj<nLevel2; rowj++) {
            tempLabels.push_back(this->rownames[rowi]+"%"+K2->rownames[rowj]);
         }
      }
      // Swap the old data with the new data.
      // Note the old data is removed when tempEvecs, tempEvals and tempLabels here go out of scope.
      this->swap(tempEvecs);
      weights->swap(tempEvals);
      std::swap(rownames,tempLabels);
   }

   std::vector<std::string> rownames, colnames;
   simpleDblVector weights;
   
};

#endif /* kernelMatrix_h */
