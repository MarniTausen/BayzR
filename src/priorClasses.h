//
//  GenericPrior.hpp
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef GenericPrior_h
#define GenericPrior_h

#include <Rcpp.h>
#include "parseColNames.h"
#include "rbayzExceptions.h"

/********

 GenericPrior object holds 'generic' prior information from R bayzPrior object, and
 has methods to sample various types of parameters: samplevar (variance), samplefreq (frequency), ...
 There is now only one constructor that needs the model-term data-column, and extracts the
 prior information from it (if available).
 Every model-term has a GenericPrior member, and the modelTerm constructor also constructs
 GenericPrior with the data-column input.
 It is not necessary that a data-column has a prior, this can be OK if there is no prior needed
 or when a default can be used.
 Internal storage is quite similar as in the R bayzPrior object, with a list of the parameters
 and other settings of the prior, and a "dist" member (from the "dist" attribute in bayzPrior)
 that has a little text-string with the name of the distribution.
 The GenericPrior constructor checks if there is a "prior" attribute in the data-column, if
 it is of type bayzPrior and list, and if it has a "dist" attribute, and if all OK copies
 the list in 'members' and dist in 'dist'.
 A status is set in the GenericPrior:
   0 (OK)
   1: there was no prior (can be OK)
   2: there was a prior attribute in the data column, but it was not an R class bayzPrior and list
   3: the prior attribute did not have a "dist" attribute (wrongly built).
 Status >1 should be reason to terminate running, but I need to think how to report errors and quit.
 Throw? 
 ********/

class GenericPrior {
   
public:
   
   GenericPrior(Rcpp::DataFrame &d, size_t c) {
      status=0;
      dist="";
      modelName = getWrapName(d, c) + "." + getVarName(d, c);
      Rcpp::RObject col = d[c];  // data column referred to as RObject
      if(col.hasAttribute("prior")) {
         Rcpp::RObject temp = col.attr("prior");
         if (temp.inherits("list") && temp.inherits("bayzPrior")) {
            if(temp.hasAttribute("dist")) {
               dist = Rcpp::as<std::string>(temp.attr("dist"));
               elements = Rcpp::as<Rcpp::List>(temp);
            }
            else {          // prior was missing dist attribute
               throw(generalRbayzError(std::string("Prior is not a correct bayzPrior object in model-term ")+modelName));
               status=3;
            }
         }
         else {            // prior attribute not of correct R class
            throw(generalRbayzError(std::string("Prior is not a correct bayzPrior object in model-term ")+modelName));
            status=2;
         }
      }
      else {               // no "prior" attribute in the col
         status=1;
      }
   }
   
   ~GenericPrior() {}
   int status;
   std::string dist, modelName;
   Rcpp::List elements;
   
   double samplevar(double ssq, size_t n) {
      if (status==1 || dist=="ichi") {
         double scale, dfprior, dftotal;
         if (status==1) {  // using uniform if no prior given
            scale = 0.0;
            dfprior = -2.0;
         }
         else {           // here not sure what will happen if "scale" or "df" are not in the list ...
            scale = elements["scale"];
            dfprior = elements["df"];
         }
         ssq += dfprior*scale;
         dftotal = dfprior+double(n);
         return(ssq/R::rchisq(dftotal));
      }
      else {
         std::string s = "Cannot sample variance with prior distribution <" + dist + "> in model-term " + modelName;
         throw(generalRbayzError(s));
      }
   }


};

#endif /* GenericPrior_h */
