//
//  data_kernel.hpp
//  Data storage class to store a 'kernel' or similarity/relationship matrix.
//  These kernels are now always stored as sets of eigenvectors and are imported
//  via ranf or ran2f functions with arguments V=. The ranf and ran2f already
//  prepare the eigenvector transform.
//  This is an implementation of dataMatrix class filling 'data' with the
//  set of eigen-vectors of the kernel, and adding extra variables eval, nEvalUsed, rrankpct.
//
//  Created by Luc Janss on 05/03/2020.
//

#ifndef dataKernel_h
#define dataKernel_h

#include <stdio.h>
#include <Rcpp.h>
#include "dataMatrix.h"

class dataKernel : public dataMatrix {

public:

   dataKernel(Rcpp::RObject col) : dataMatrix(col) {

      // The exceptions here should not occur ... as long as set-up of ranf and ran2f and code
      // in rbayz main do not change. It may be possible to add name of the column in error msgs?
      if (col.hasAttribute("evectors"))
         data = Rcpp::as<Rcpp::NumericMatrix>(col.attr("evectors"));
      else
         throw(generalRbayzError("Eigen-decomp storage is corrupted"));
      if (col.hasAttribute("evalues"))
         eval = Rcpp::as<Rcpp::NumericVector>(col.attr("evalues"));
      else
         throw(generalRbayzError("Eigen-decomp storage is corrupted"));
      if(data.nrow() != data.ncol())
         throw(generalRbayzError("Eigen-decomp has non-square e-vectors matrix"));
      if(eval.size() != data.ncol())
         throw(generalRbayzError("Eigen-decomp vectors and values size do not match"));
      if(col.hasAttribute("rrankpct"))
         rrankpct = col.attr("rrankpct");
      else
         rrankpct = 99;
      double sumeval = 0.0l;
      for (size_t i = 0; i<eval.size() && eval[i] > 0; i++) sumeval += eval[i];
      double eval_cutoff = rrankpct * sumeval / 100.0l;
      nEvalUsed = 0;
      sumeval = 0.0l;
      while (sumeval < eval_cutoff) sumeval += eval[nEvalUsed++];
      // this message can be moved to the modelTerm routine; and it must be differen for ran2f_2cor.
      Rcpp::Rcout << "In ranf with V rrankpct=" << rrankpct << " uses " << nEvalUsed << "eigenvectors\n";

   /* version to use low-level C arrays
      Rcpp::NumericMatrix evecR  // local variable that only lives in the constructor
            = Rcpp::as<Rcpp::NumericMatrix>(col.attr("evectors"));
       ... find out sizes to use (can already check evals and only copy the needed ones)
      evec = static_cast<double *> (malloc(N*sizeof(double)));
      if (evec==NULL) throw(generalRbayzError("Out of memory"));
      ... copy data from evecR in evec
    */
   }

   ~dataKernel() {
       // needs clean up when using new or malloc in constructor
   }
   
   void check(){}

   size_t nEvalUsed;
   double rrankpct;
};

#endif /* dataKernel_h */
