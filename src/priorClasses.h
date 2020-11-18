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
#include "parseFunctions.h"
#include "rbayzExceptions.h"

/* GenericPrior object: to hold prior information in a generic way, bascically
   a map of (parameter, value) pairs and a "dist" string for the distribution.
   However, the construction has moved from R to inside bayz, so that
   it needs a check that the `dist` is compatible with the specified parameters.
   Better make separate classes for each prior? But how to initialise the correct object?
 */

class GenericPrior {
   
public:
   
   // initialisation from a model term, a string like "rf(A:B, V=..., prior=ichi(...))"
   GenericPrior(std::string modelTerm) {
      size_t pos1, pos2, pos3;
      if ( (pos1=modelTerm.find("prior=")) != std::string::npos) {
         // Note: if prior= is not found, nothing happens here and a default will be used.
         useDefault=FALSE;
         pos1 += 6;   // now pos1 has moved where there should be the distribution name
         pos2 = modelTerm.find('(',pos1);
         pos3 = modelTerm.find(')',pos2);
         dist = modelTerm.substr(pos1,(pos2-pos1));
         std::vector<std::string> parlist = splitString(modelTerm.substr(pos2+1,(pos3-pos2-1)),',');
         std::string parname;
         double parvalue;
         for(size_t i=0; i<parlist.size(); i++) {
            if((pos1 = parlist[i].find('=')) == std::string::npos) {
               throw(generalRbayzError("Error in parameter specification <"
                                       +parlist[i]+"> in "+modelTerm));
            }
            else {
               parname = parlist[i].substr(0,pos1);
               try {
                  parvalue = std::stod(parlist[i].substr(pos1+1,std::string::npos));
               }
               catch (std::exception& e) {
                  throw(generalRbayzError("Error to get value from parameter <"
                                          +parlist[i]+"> in "+modelTerm));
               }
               param.insert(std::make_pair(parname,parvalue));
            }
         }
         if(dist=="ichi") {
            if(param.find("scale") == param.end())
               throw(generalRbayzError("Error ichi prior is missing scale parameter in "+modelTerm));
            if(param.find("df") == param.end())
               throw(generalRbayzError("Error ichi prior is missing df parameter in "+modelTerm));
         }
         // else if(dist==...)  // add checks for other distributions
         else {
            throw(generalRbayzError("Unknown prior distribution in "+modelTerm));
         }
      }
   }

   ~GenericPrior() {}

   bool useDefault=TRUE;
   std::string dist, modelName;
   std::map<std::string,double> param;
   
   double samplevar(double ssq, size_t n) {
      if (useDefault || dist=="ichi") {
         double scale, dfprior, dftotal;
         if (useDefault) {  // using uniform if no prior given
            scale = 0.0;
            dfprior = -2.0;
         }
         else {           // here not sure what will happen if "scale" or "df" are not in the list ...
            scale = param.find("scale")->second;
            dfprior = param.find("df")->second;
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
