//
//  dataFactor.h
//  rbayz
//
//  Created by Luc Janss on 03/08/2018.
//

#ifndef dataFactor_h
#define dataFactor_h

#include <Rcpp.h>
#include <vector>
#include <string>

void CharVec2cpp(std::vector<std::string> & labels, Rcpp::CharacterVector templabels);

class dataFactor {
   
public:

   // It is now possible to setup an 'empty' dataFactor object, it can be
   // filled with the addVariable() function.
   dataFactor() {
      Nvar=0;
   }
   
   // Constructors to setup one or first factor, with three interfaces, but
   // in the end all use a function that handles an RObject.
   // The string-interface is interpreted as passing the name of an object
   // that should be in the R environment, and names are already checked in
   // the modelBase constructor.
   dataFactor(Rcpp::DataFrame &d, size_t col) {
      Rcpp::RObject Rcol = Rcpp::as<Rcpp::RObject>(d[col]);
      setupFirstVariable(Rcol);
   }

   dataFactor(std::string varname) {
      Rcpp::Environment Renv;
      Rcpp::RObject Rcol = Renv[varname];
      setupFirstVariable(Rcol);
   }

   dataFactor(Rcpp::RObject Rcol) {
      setupFirstVariable(Rcol);
   }

   void setupFirstVariable(Rcpp::RObject col) {
      if (!Rf_isFactor(col)) {
         throw generalRbayzError("Variable is not a factor (unfortunately cannot get the name here)\n");
      }
      data = Rcpp::as<Rcpp::IntegerVector>(col);
      for(size_t i=0; i < data.size(); i++) {
         Rcpp::Rcout << " " << data[i];
      }
      for (size_t i=0; i<data.size(); i++)
         data[i] -= 1;
      Rcpp::Rcout << "\n";
      Rcpp::CharacterVector templabels = data.attr("levels");
      CharVec2cpp(labels, templabels);
      Nvar=1;
      Rcpp::Rcout << "Factor set-up with " << labels.size() << " levels\n";
      for(size_t i=0; i < data.size(); i++) {
         Rcpp::Rcout << " " << data[i];
      }
      Rcpp::Rcout << "\n";
      for(size_t i=0; i < labels.size(); i++) {
         Rcpp::Rcout << " " << labels[i];
      }
      Rcpp::Rcout << "\n";

   }

   // addVariables adds another variable in a factor
   void addVariable(Rcpp::DataFrame &d, size_t col) {
      Rcpp::RObject Rcol = Rcpp::as<Rcpp::RObject>(d[col]);
      addVariable(Rcol);
   }

   void addVariable(std::string varname) {
      Rcpp::Environment Renv;
      Rcpp::RObject Rcol = Renv[varname];
      addVariable(Rcol);
   }

   void addVariable(Rcpp::RObject Rcol) {
      if(Nvar==0)
         setupFirstVariable(Rcol);
      else { // add (interact) another factor with already stored factor(s)
         if (!Rf_isFactor(Rcol)) {
            throw generalRbayzError("Variable is not a factor (unfortunately cannot get the name here)\n");
         }
         std::vector<std::string> oldlabels(labels);
         std::vector<std::string> newFactorLabels;
         CharVec2cpp(newFactorLabels, data.attr("levels"));
         Rcpp::IntegerVector newdata = Rcpp::as<Rcpp::IntegerVector>(Rcol);
         for (size_t i=0; i<newdata.size(); i++)
            newdata[i] -= 1;
         size_t nLevel1=oldlabels.size();
         size_t nLevel2=newFactorLabels.size();
         labels.resize(nLevel1 * nLevel2);
         for(size_t i=0; i<nLevel1; i++) {  // generate the new labels
            for(size_t j=0; j<nLevel2; j++) {
               labels[i*nLevel2+j] = oldlabels[i] + "%" + newFactorLabels[j];
            }
         }
         // Replace existing data-level-coding with codes to match the new interaction
         for(size_t i=0; i<data.size(); i++) {
            data[i] = data[i]*nLevel2 + newdata[i];
         }
         Nvar++;
      }
   }

   ~dataFactor() {
   }

   Rcpp::IntegerVector data;
   std::vector<std::string> labels;
   int Nvar;  // The number of (interacting) variables in this factor

};

#endif /* dataFactor_h */
