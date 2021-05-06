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
//#include "nameTools.h"  // strange, including nameTools.h does not work to make
                          // addMatrixNames available here ...

int addMatrixNames(std::vector<std::string> & names, Rcpp::NumericMatrix & mat, int dim);

class kernelMatrix : public simpleMatrix {

public:

   kernelMatrix(Rcpp::RObject col, std::string & name) : simpleMatrix(), weights() {
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
      double rrankpct;
      size_t nColUsed = 0;
      if(col.hasAttribute("rrankpct"))
         rrankpct = col.attr("rrankpct");
      else
         rrankpct = 99;
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
      addMatrixNames(rownames, kerneldata, 1);
      // if want to add colnames, it can be evec1, evec2, etc. The original colnames
      // on the kernel (if present) have become meaningless after the eigen-decomp.
   }

   ~kernelMatrix() {
   }

   std::vector<std::string> rownames, colnames;
   simpleDblVector weights;
   
};

#endif /* kernelMatrix_h */
