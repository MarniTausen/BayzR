#include <vector>
#include <string>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::plugins("cpp11")]]

#include "parseColNames.h"
#include "modelBase.h"
#include "modelResp.h"
#include "modelMean.h"
#include "modelFixf.h"
#include "modelRanf.h"
#include "modelFreg.h"
#include "modelRanf_cor.h"
#include "modelRan2f.h"
#include "modelRan2f_2cor.h"
#include "rbayzExceptions.h"

// These functions are defined below the main function
std::vector<std::string> modelTerms = parseModel(Rcpp::String modelFormula);
void buildModelTerm(Rcpp::DataFrame & modelFrame, size_t col, std::vector<modelBase *> & model, Rcpp::RObject &terms);
void collectParInfo(std::vector<modelBase *> & model, Rcpp::CharacterVector & parNames,
                    Rcpp::LogicalVector & parHyper, Rcpp::IntegerVector & parSizes,
                    Rcpp::IntegerVector & parEstFirst, Rcpp::IntegerVector & parEstLast,
                    Rcpp::IntegerVector & parModelNr, Rcpp::IntegerVector & parLogged,
                    Rcpp::CharacterVector & parLoggedNames, Rcpp::CharacterVector & estimNames);
void writeLoggedSamples(size_t & cycle, std::vector<modelBase *> & model, Rcpp::IntegerVector & parLogged,
                        Rcpp::CharacterVector & parLoggedNames, Rcpp::IntegerVector & parModelNr, bool silent);
void collectPostStats(std::vector<modelBase *> & model, Rcpp::NumericVector & postMean,
                      Rcpp::NumericVector & postSD);
void collectLoggedSamples(std::vector<modelBase *> &model, Rcpp::IntegerVector & parModelNr,
                          Rcpp::IntegerVector & parLogged, Rcpp::NumericMatrix & loggedSamples, size_t save);

// The interface defines to pass modelFrame by value, this is a quite 'light' copy (a bunch
// of pointers); also modelFrame is modified within the main function, but these modifications
// don't need to persist. For sub-functions, the modelFrame is passed around by reference.
// The return value is a List, to check if program terminated normally or with errors check $nError.

// [[Rcpp::export]]
Rcpp::List rbayz_cpp(Rcpp::String modelFormula, Rcpp::DataFrame modelFrame,
                     Rcpp::IntegerVector chain, bool silent) {

   // Some check of chain settings is needed. Also the rbayz wrapper function now
   // handles chain being NULL, but it can be done here, so that warning message
   // can come in output.
   
   // Vectors to collect error and 'notes' messages are defined outside the try-block,
   // so it remains available in case of errors.
   Rcpp::CharacterVector errorMessages;
   Rcpp::CharacterVector notesMessages;
   std::string lastDone;
   std::vector<std::string> modelTerms = parseModel(modelFormula);

   try {     // a large try-block wraps nearly all of code, in case of normal exectution
             // the code builds a return list and returns before the catch().
             // In case of errors, catch() builds a return list with the messages vector.

      // I found a strange behaviour with modelFrame (a DataFrame): as soon as I add a column to
      // it (makes a copy), R converts it to a list, and attributes row.names, and terms get lost,
      // and method .nrow() is no longer available! So I need to collect these things before losing them.
      size_t nRow = modelFrame.nrow();
      size_t nCol = modelFrame.size();                  // nr columns of input, later it will grow
      Rcpp::RObject terms = modelFrame.attr("terms");
      Rcpp::CharacterVector rowNames = modelFrame.attr("row.names");
      if(terms.isNULL()) throw(generalRbayzError("Input model frame is lacking terms"));
      
      // All model-terms must have access to the vector of residuals and residual precisions.
      // To organize this, two vectors are added to the model-frame.
      Rcpp::NumericVector e(nRow,0);
      Rcpp::NumericVector v(nRow,1);
      modelFrame.push_back(e,"residual");
      modelFrame.push_back(v,"residPrec");

      // A vector of pointers to modelBase objects
      std::vector<modelBase *> model;

      // Build model by building a model-term from each input column (nCol columns) in the model-frame.
      // Note: when model-frame has flagged intercept, a 'mean' term is added when handling col zero.
      bool modelBuildError=FALSE;
      for(size_t col=0; col<nCol; col++) {
         try {
            buildModelTerm(modelFrame, col, model, terms);
         }
         catch (generalRbayzError &err) {        // model building continues after a generalRBayzError,
            errorMessages.push_back(err.what());      // modelBuildError is set TRUE and will throw exception below.
            modelBuildError=TRUE;                // Any other exception will directly jump to outer catch().
         }
      }
      if (modelBuildError)
         throw(generalRbayzError("Bayz terminates after model building"));

      // Parameter information vectors
      Rcpp::CharacterVector parNames;               // collected names of used hpar and par vectors (where size>0)
      Rcpp::LogicalVector parHyper;                 // if it is a hpar (TRUE) or par (FALSE) vector
      Rcpp::IntegerVector parSizes;                 // size (number of levels) in each hpar and par vector
      Rcpp::IntegerVector parEstFirst,parEstLast;   // First and last number in estimates vector (R coding from 1)
      Rcpp::IntegerVector parModelNr;               // which model-number the parameter vector comes from
      Rcpp::IntegerVector parLogged;                // If the parameter will be logged
      Rcpp::CharacterVector parLoggedNames;         // Shortened list of names of logged parameters
      Rcpp::CharacterVector estimNames;             // Names of all 'estimates' (individual levels)
      
      collectParInfo(model, parNames, parHyper, parSizes, parEstFirst, parEstLast, parModelNr,
                     parLogged, parLoggedNames, estimNames);
      size_t nEstimates = estimNames.size();

      // Make list of output sample cycle-numbers; this also checks if the chain settings actually make 1 or more output.
      Rcpp::IntegerVector outputCycleNumbers;
      for (size_t cycle=1; cycle <= chain[0]; cycle++) {      // exactly the same loop and condition to make
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0))  // output samples below
            outputCycleNumbers.push_back(cycle);
      }
      size_t nSamples = outputCycleNumbers.size();
      if (nSamples==0) throw (generalRbayzError("The chain settings do not make any output"));
      Rcpp::NumericMatrix loggedSamples(int(nSamples),int(parLoggedNames.size()));
      Rcpp::NumericVector postMean(nEstimates,0);
      Rcpp::NumericVector postSD(nEstimates,0);
      Rcpp::colnames(loggedSamples) = parLoggedNames;
      int nShow = chain[0]/5;
      lastDone="Preparing to run MCMC";

      // Run the model by calling the sample() method for each modelTerm
      for (size_t cycle=1, save=0; cycle <= chain[0]; cycle++) {
         for(size_t mt=0; mt<model.size(); mt++) model[mt]->sample();
         if (cycle % nShow == 0 )
            writeLoggedSamples(cycle, model, parLogged, parLoggedNames, parModelNr, silent);
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0) ) {
            for(size_t mt=0; mt<model.size(); mt++) model[mt]->prepForOutput();
            collectPostStats(model, postMean, postSD);
            collectLoggedSamples(model, parModelNr, parLogged, loggedSamples, save);
            save++;  // save is counter for output (saved) cycles
         }
      }
      lastDone="Finished running MCMC";

      // Build result list for normal termination
      Rcpp::List result = Rcpp::List::create();
      Rcpp::DataFrame parInfo = Rcpp::DataFrame::create
               (Rcpp::Named("ModelNr")=parModelNr, Rcpp::Named("Hyper")=parHyper,
                Rcpp::Named("Size")=parSizes, Rcpp::Named("EstStart")=parEstFirst,
                Rcpp::Named("EstEnd")=parEstLast, Rcpp::Named("Logged")=parLogged);
      parInfo.attr("row.names") = parNames;
      lastDone="Setting up return list";
      for(size_t i=0; i<nEstimates; i++)
          postMean[i] /= double(nSamples);   // finish computing post mean and SD, so far it is sum and sum squares
      for(size_t i=0; i<nEstimates; i++)
          postSD[i] = sqrt(postSD[i]/double(nSamples) - postMean[i]*postMean[i]);
      lastDone="Computing postMeans and PostSDs";
      Rcpp::DataFrame estimates = Rcpp::DataFrame::create
               (Rcpp::Named("postMean")=postMean, Rcpp::Named("postSD")=postSD);
      estimates.attr("row.names") = estimNames;
      lastDone="Setting up estimates dataframe";
      result.push_back(0,"nError");
      result.push_back(parInfo,"Parameters");
      result.push_back(loggedSamples,"Samples");
      result.push_back(estimates,"Estimates");
      lastDone="Filling return list";
      return(result);   // normal termination

   } catch (std::exception &err) {  // this includes generalRbayzErrors
      errorMessages.push_back(err.what());
   }
   catch (...) {
      errorMessages.push_back("An unknown error occured in bayz after: "+lastDone);
   }

   // Note: program flow only comes here in case of errors (normal return is above);
   // build a return list that only has the error messages list.
   Rcpp::List result = Rcpp::List::create();
   result.push_back(errorMessages.size(),"nError");
   result.push_back(errorMessages,"Errors");
   Rcpp::Rcout << "Bayz finished with errors - use summary() or check $Errors" << std::endl;
   return(result);

}

void buildModelTerm(Rcpp::DataFrame & modelFrame, size_t col, std::vector<modelBase *> & model, Rcpp::RObject &terms) {
   std::string s = getWrapName(modelFrame, col);
   if (s=="" && col==0) {
      model.push_back(new modelResp(modelFrame, col));
      int k=0;
      if (terms.hasAttribute("intercept"))
         k = terms.attr("intercept");
      if (k==1)   // model has intercept
         model.push_back(new modelMean(modelFrame, 0));
   }
   else if (s=="fixf")
      model.push_back(new modelFixf(modelFrame, col));
   else if(s=="ranf") {
      Rcpp::RObject thiscol = modelFrame[col];
      if(thiscol.hasAttribute("evalues")) {
         model.push_back(new modelRanf_cor(modelFrame, col));
      }
      else {
         model.push_back(new modelRanf(modelFrame, col));
      }
   }
   else if(s=="ran2f") {
      Rcpp::RObject thiscol = modelFrame[col];
      Rcpp::RObject secondcol = thiscol.attr("factor2");
      if(thiscol.hasAttribute("evectors") && secondcol.hasAttribute("evectors")) {
         model.push_back(new modelRan2f_2cor(modelFrame, col));
      }
      else if (!thiscol.hasAttribute("evectors") && !secondcol.hasAttribute("evectors")) {
         model.push_back(new modelRan2f(modelFrame, col));
      }
      else if (thiscol.hasAttribute("evectors") && !secondcol.hasAttribute("evectors")) {
         throw(generalRbayzError("Version of ran2f not yet implemented"));
      }
      else {
         throw(generalRbayzError("In ran2f(...,V1=,V2=) cannot set V2 only, swap the variables and set V1"));
      }
   }
   else if (s=="freg")
      model.push_back(new modelFreg(modelFrame, col));
   else if (col==0 && s!="") {
      throw(generalRbayzError("Cannot handle wrapper function on response column"));
   }
   else {
      throw(generalRbayzError("Unknown wrapper function on data column"));
   }
   return;
}

void collectParInfo(std::vector<modelBase *> & model, Rcpp::CharacterVector & parNames,
                    Rcpp::LogicalVector & parHyper, Rcpp::IntegerVector & parSizes,
                    Rcpp::IntegerVector & parEstFirst, Rcpp::IntegerVector & parEstLast,
                    Rcpp::IntegerVector & parModelNr, Rcpp::IntegerVector & parLogged,
                    Rcpp::CharacterVector & parLoggedNames, Rcpp::CharacterVector & estimNames) {

   // First collect list of used hpar and par vectors, and set what model(-term) they belong to.
   for(size_t mt=0; mt<model.size(); mt++) {     // Done in two loops over model-terms, all hpar will come
      if ( model[mt]->hpar.size() > 0 ) {        // first, then all par vectors
         parHyper.push_back(TRUE);
         parModelNr.push_back(mt);
         parSizes.push_back(model[mt]->hpar.size());
         parNames.push_back(model[mt]->hparName);
         if (model[mt]->hpar.size()==1) {        // log all hyper-parameters with 1 level
            parLogged.push_back(1);
            parLoggedNames.push_back(model[mt]->hparName);
         }
         else
            parLogged.push_back(0);
      }
   }
   for(size_t mt=0; mt<model.size(); mt++) {
      if (model[mt]->par.size() > 0 ) {
         parHyper.push_back(FALSE);
         parModelNr.push_back(mt);
         parSizes.push_back(model[mt]->par.size());
         parNames.push_back(model[mt]->parName);
         if (model[mt]->par.size()==1) {       // log all parameters with 1 level. Could also consider to log
            parLogged.push_back(2);            // 2nd level of 2-level fixf parameters (code 3 ) ...
            parLoggedNames.push_back(model[mt]->parName);
         }
         else
            parLogged.push_back(0);
      }
   }
   
   // Fill names of all estimates (parameter%level), and fill parEstFirst, parEstLast which tell where each
   // parameter-vector can be found in the estimates list. This runs one loop, but that requires to set pointers to
   // par or hpar vectors (names, levels); alternative is to copy 2nd part of code two times in loops above....
   size_t nEstimates=1;  // to fill First and Last on index-base 1, after loop nEstimates
                         // will be 1 more than total levels
   {
      std::vector<double> *par_ptr;
      std::string *name_ptr;
      std::vector<std::string> *label_ptr;
      std::string s;
      for(size_t par=0; par<parModelNr.size(); par++) {
         if (parHyper[par]) {
            par_ptr = &(model[parModelNr[par]]->hpar);
            name_ptr = &(model[parModelNr[par]]->hparName);
            label_ptr = &(model[parModelNr[par]]->hparLabels);
         }
         else {
            par_ptr = &(model[parModelNr[par]]->par);
            name_ptr = &(model[parModelNr[par]]->parName);
            label_ptr = &(model[parModelNr[par]]->parLabels);
         }
         parEstFirst.push_back(nEstimates);
         nEstimates += par_ptr->size();
         parEstLast.push_back(nEstimates-1);
         if (par_ptr->size()==1) {  // for parameter-vectors size 1, estimate name is same as parameter name
            estimNames.push_back(*name_ptr);
         }
         else {                    // parameter-vector >1: the estimNames get appended %level
            if (par_ptr->size() != label_ptr->size()) {
               // Here is a piece of code filling dummy labels (1,2,...) when the length of
               // the labels does not match the length of the parameter vector. This mostly happens
               // when label-vector is not filled, but also avoids conflict when the two don't match.
               for (size_t level=0; level < par_ptr->size(); level++) {
                  s = *name_ptr + "%" + std::to_string(level+1);
                  estimNames.push_back(s);
               }
            }
            else {
               for(size_t level=0; level < label_ptr->size(); level++) {
                  s = *name_ptr + "%" + (*label_ptr)[level];
                  estimNames.push_back(s);
               }
            }
         }
      }
   }
   
   if (estimNames.size() != (nEstimates-1))
      throw(generalRbayzError("Error in collecting parameter-info lists"));
      // This error should in principle not happen, every time nEstimates is increased with the
      // size of a parameter vector, also this same number of names is pushed_back in estimNames.

   return;

}

// This is the only part that writes to the screen under normal operation (if there are no errors)
void writeLoggedSamples(size_t & cycle, std::vector<modelBase *> & model, Rcpp::IntegerVector & parLogged,
                        Rcpp::CharacterVector & parLoggedNames, Rcpp::IntegerVector & parModelNr, bool silent) {
   if (silent) return;
   static int writtenHead=0;
   double x=0.0;
   if (!writtenHead) {
      Rcpp::Rcout << "cycle ";
      for(size_t i=0; i<parLoggedNames.size(); i++)
      Rcpp::Rcout << parLoggedNames[i] << " ";
      Rcpp::Rcout <<  std::endl;
      writtenHead=1;
   }
   Rcpp::Rcout << cycle << " ";
   for(size_t i=0; i<parLogged.size(); i++) {
      if (parLogged[i] > 0) {
         if (parLogged[i] == 1) x = model[parModelNr[i]]->hpar[0];
         else if (parLogged[i] == 2) x = model[parModelNr[i]]->par[0];
         else if (parLogged[i] == 3) x = model[parModelNr[i]]->par[1];
         Rcpp::Rcout << x << " ";
      }
   }
   Rcpp::Rcout << std::endl;
}

void collectPostStats(std::vector<modelBase *> & model, Rcpp::NumericVector & postMean,
                      Rcpp::NumericVector & postSD) {
   size_t k=0;
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->hpar.size() > 0 ) {
         for(size_t i=0; i<model[mt]->hpar.size(); i++) {
            postMean[k] += model[mt]->hpar[i];
            postSD[k] += (model[mt]->hpar[i])*(model[mt]->hpar[i]);
            k++;
         }
      }
   }
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->par.size() > 0 ) {
         for(size_t i=0; i<model[mt]->par.size(); i++) {
            postMean[k] += model[mt]->par[i];
            postSD[k] += (model[mt]->par[i])*(model[mt]->par[i]);
            k++;
         }
      }
   }
}

void collectLoggedSamples(std::vector<modelBase *> &model, Rcpp::IntegerVector & parModelNr,
                          Rcpp::IntegerVector & parLogged, Rcpp::NumericMatrix & loggedSamples, size_t save) {
   size_t k=0;
   double x=0.0;
   if(save > (loggedSamples.nrow()-1)) {
      throw(generalRbayzError("Storage for logged samples out of bounds"));
   }
   for(size_t i=0; i<parLogged.size(); i++) {
      if (parLogged[i] > 0) {
         if (parLogged[i] == 1) x = model[parModelNr[i]]->hpar[0];
         else if (parLogged[i] == 2) x = model[parModelNr[i]]->par[0];
         else if (parLogged[i] == 3) x = model[parModelNr[i]]->par[1];
         loggedSamples(save,k) = x;
         k++;
      }
   }
}
