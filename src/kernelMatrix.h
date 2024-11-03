//  kernelMatrix.h
//  Storage of kernel as eigen-decomposition.
//  Derives from simpleMatrix using the simpleMatrix() empty contructor, so that
//  no matrix data is stored yet from the parent constructor.
//
//  Created by Luc Janss on 06/05/2021.
//

#ifndef kernelMatrix_h
#define kernelMatrix_h

#include <stdio.h>
#include <Rcpp.h>
#include "rbayzExceptions.h"
#include "parseFunctions.h"
#include <string>
#include <vector>
#include "labeledMatrix.h"
#include "simpleVector.h"

extern std::vector<std::string> Messages;
extern bool needStop;

class kernelMatrix : public labeledMatrix {

public:

   kernelMatrix(Rcpp::RObject col, std::string & name,    std::map<std::string, std::string> & mtoptions)
         : labeledMatrix(), weights() {
//    note: kernelMatrix starts with an empty labeledMatrix, parent constructors have not done anything,
//    accept for having the matrix and vectors for storing data and row and column labels.
//    The initWith() used at the end to copy eigenvec contents in the object is not simpleMatrix' initWith,
//    but labeledMatrix' version that also handles copying row and column labels.
      Rcpp::NumericMatrix kerneldata = Rcpp::as<Rcpp::NumericMatrix>(col);
   	Rcpp::Function eig("eigen");
	   Rcpp::List eigdecomp;
   	try {
	   	eigdecomp = eig(kerneldata);
   	}
   	catch(std::exception &err) {
         throw(generalRbayzError("An error occurred running eigen(): "+std::string(err.what())));
   	}
   	Rcpp::NumericVector eigvalues = eigdecomp["values"];
	   Rcpp::NumericMatrix eigvectors = eigdecomp["vectors"];
      // re-attach the dimnames again to the eigvectors matrix for correct further processing
      if (kerneldata.hasAttribute("dimnames")) {
         Rcpp::List dimnames = Rcpp::as<Rcpp::List>(kerneldata.attr("dimnames"));
         eigvectors.attr("dimnames") = dimnames;
      }
      // Determine how many evectors to keep
      std::string dimopt, dimpopt;
      double dim_pct=0;
      int dim_size=0;
      if( ( dimopt = mtoptions["dim"]) != "") {
         dim_size = str2int(dimopt, ("dim="+dimopt));
         if(dim_size <= 0 || dim_size > eigvalues.size()) {
            Messages.push_back("Warning: invalid dim setting <"+dimopt+"> processing kernel "+name+", setting default dimp=80");
            dim_size=0;  // if dim not well set this does not trigger error,
            dim_pct=80;  // but goes back to cutting off on 80% of variance.
         }
      }
      else {  // without 'dim' option, check for 'dimp' (note: dim will be used when both are set!)
         if( (dimpopt = mtoptions["dimp"]) != "") {
            dim_pct = str2dbl(dimpopt,("dimp="+dimpopt));
            if(dim_pct <= 0 || dim_pct > 100) {
               Messages.push_back("Warning: invalid dimp setting <"+dimpopt+"> processing kernel "+name+", setting default dimp=80");
               dim_pct=80;  // also here not error, but fo back to default dimp=80
            }
         }
         else  // no options set: take default
            dim_pct=80;
      }
      if(dim_size==0) {      // cut-off on %variance, need to find the matching number of evecs
         double sumeval = 0.0l;                             //  \/ only summing the positive ones!
         for (size_t i = 0; i < unsigned(eigvalues.size()) && eigvalues[i] > 0; i++)
            sumeval += eigvalues[i];
         double eval_cutoff = dim_pct * sumeval / 100.0l;
         sumeval = 0.0l;
         while (sumeval < eval_cutoff) sumeval += eigvalues[dim_size++];
         std::string s = "Note: for kernel " + name + " using dimp=" + std::to_string(dim_pct) + " takes "
                      + std::to_string(dim_size) + " eigenvectors";
         Messages.push_back(s);
      }
      this->initWith(eigvectors, name, dim_size);
      weights.initWith(eigvalues, dim_size);
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
               tempEvals.data[k] = this->weights[i] * K2->weights[j];
               for(size_t rowi=0; rowi<nLevel1; rowi++) {
                  for(size_t rowj=0; rowj<nLevel2; rowj++) {
                     tempEvecs.data[k][i*nLevel2+j] = this->data[i][rowi] * K2->data[j][rowj];
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
            tempLabels.push_back(this->rownames[rowi]+"."+K2->rownames[rowj]);
         }
      }
      Rcpp::Rcout << "Interaction kernel retains " << nEvalUsed << " eigenvalues\n";
      // Swap the old data with the new data.
      // Note the old data is removed when tempEvecs, tempEvals and tempLabels here go out of scope.
      this->swap(&tempEvecs);
      weights.swap(&tempEvals);
      std::swap(rownames,tempLabels);
   }

   simpleDblVector weights;
   
};

#endif /* kernelMatrix_h */
