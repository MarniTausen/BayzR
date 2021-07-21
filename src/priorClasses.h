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

/* GenericPrior object: to hold prior information in a generic way, basically
   a map of (parameter, value) pairs and a "dist" string for the distribution.
   However, the construction has moved from R to inside bayz, so that
   it needs a check that the `dist` is compatible with the specified parameters.
   Better make separate classes for each prior? But how to initialise the correct object?
 */

class GenericPrior {
   
public:
   
   // initialisation from the string after prior=, or empty
   GenericPrior(std::string priordescr) {
      size_t pos1, pos2, pos3;
      if ( priordescr != "") {
         // Note: if priordescr is empty nothing happens here and a default will be used.
         useDefault=FALSE;
         pos2 = priordescr.find('(');         // syntax is like ichi(par1=val1, par=val2, etc),
         pos3 = priordescr.find(')',pos2);    // here want to get the part before first parenthesis
         dist = priordescr.substr(0,pos2);    // as the dist name, and split part inside parenthesis
         std::vector<std::string> parlist     // in a list of par=value strings
                 = splitString(priordescr.substr(pos2+1,(pos3-pos2-1)),",");
         std::string parname;
         double parvalue;
         for(size_t i=0; i<parlist.size(); i++) {
            if((pos1 = parlist[i].find('=')) == std::string::npos) {
               throw(generalRbayzError("Error in parameter specification <"
                                       +parlist[i]+"> in "+priordescr));
            }
            else {
               parname = parlist[i].substr(0,pos1);
               try {
                  parvalue = std::stod(parlist[i].substr(pos1+1,std::string::npos));
               }
               catch (std::exception& e) {
                  throw(generalRbayzError("Error to get value from parameter <"
                                          +parlist[i]+"> in "+priordescr));
               }
               param.insert(std::make_pair(parname,parvalue));
            }
         }
         if(dist=="ichi") {
            if(param.find("scale") == param.end())
               throw(generalRbayzError("Error ichi prior is missing scale parameter in "+priordescr));
            if(param.find("df") == param.end())
               throw(generalRbayzError("Error ichi prior is missing df parameter in "+priordescr));
         }
         // else if(dist==...)  // add checks for other distributions
         else {
            throw(generalRbayzError("Unknown prior distribution in "+priordescr));
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
