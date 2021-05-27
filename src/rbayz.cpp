#include <vector>
#include <string>
#include <Rcpp.h>
#include <cmath>
#include "parseFunctions.h"
#include "dcModelTerm.h"
#include "modelBase.h"
#include "modelResp.h"
#include "modelMean.h"
#include "modelFixf.h"
#include "modelRanfi.h"
#include "modelFreg.h"
#include "modelRreg.h"
#include "modelRanf_cor.h"
#include "rbayzExceptions.h"
#include "simpleMatrix.h"
#include "simpleVector.h"
#include "modelVar.h"
#include "indepVarStr.h"

// [[Rcpp::plugins("cpp11")]]

// These functions are defined below the main function
void buildModel(std::vector<modelBase *> &, dcModelTerm &, modelBase*);
void buildVarModel(std::vector<modelBase *> &, dcModelTerm &, modelBase*);
void insertIndepVar(std::vector<modelBase *> &, dcModelTerm &, modelBase*);
void collectParInfo(std::vector<modelBase *> & model, Rcpp::CharacterVector & parNames,
                    Rcpp::LogicalVector & parHyper, Rcpp::IntegerVector & parSizes,
                    Rcpp::IntegerVector & parEstFirst, Rcpp::IntegerVector & parEstLast,
                    Rcpp::IntegerVector & parModelNr, Rcpp::CharacterVector & parModelFunc, Rcpp::IntegerVector & parLogged,
                    Rcpp::CharacterVector & parLoggedNames, Rcpp::CharacterVector & estimNames);
void writeLoggedSamples(size_t & cycle, std::vector<modelBase *> & model, Rcpp::IntegerVector & parLogged,
                        Rcpp::CharacterVector & parLoggedNames, Rcpp::IntegerVector & parModelNr, bool silent);
void collectPostStats(std::vector<modelBase *> & model, Rcpp::NumericVector & postMean,
                      Rcpp::NumericVector & postSD);
void collectLoggedSamples(std::vector<modelBase *> & model, Rcpp::IntegerVector & parModelNr,
                          Rcpp::IntegerVector & parLogged, Rcpp::NumericMatrix & loggedSamples, size_t save);

// [[Rcpp::export]]
Rcpp::List rbayz_cpp(Rcpp::Formula modelFormula, Rcpp::DataFrame inputData,
                     Rcpp::IntegerVector chain, int silent)
{
   // Some check of chain settings is needed. Also the rbayz wrapper function now
   // handles chain being NULL, but it can be done here, so that warning message
   // can come in output.
   
   // Vectors to collect error and 'notes' messages are defined outside the try-block,
   // so it remains available in case of errors.
   Rcpp::CharacterVector errorMessages;
   Rcpp::CharacterVector notesMessages;
   std::string lastDone;

   try {     // a large try-block wraps nearly all of code, in case of normal execution
             // the code builds a return list and returns before the catch().
             // In case of errors, catch() builds a return list with the messages vector.

      // split the modelFormula in a list of response (LHS) and explanatory (RHS) terms.
      // getModelRHSTerms makes sure there is a "0" or "1" to specify if a mean is needed.
      std::string formulaAsCppstring = convertFormula(modelFormula);
      removeSpaces(formulaAsCppstring);
      std::vector<std::string> modelLHSTerms = getModelLHSTerms(formulaAsCppstring);
      std::vector<std::string> modelRHSTerms = getModelRHSTerms(formulaAsCppstring);

      // Vectors to hold pointers to the modelling objects
      std::vector<modelBase *> model;

      // Build the modelling objects. The RHS terms are nested within the LHS (response)
      // terms to build a explanatory term for every response.
      for(size_t resp=0; resp<modelLHSTerms.size(); resp++) {
         dcModelTerm parsedResponseModelDescr(modelLHSTerms[resp], inputData);
         modelResp* tempptr = new modelResp(parsedResponseModelDescr, NULL);
         modelBase* responseModel = dynamic_cast<modelBase *>(tempptr);
         model.push_back(responseModel);
         insertIndepVar(model, parsedResponseModelDescr, model.back());
         tempptr->varModel = dynamic_cast<indepVarStr *>(model.back());
         for(size_t term=0; term<modelRHSTerms.size(); term++) {
            dcModelTerm modeldescr(modelRHSTerms[term], inputData);
            buildModel(model, modeldescr, responseModel);
            // Same for hierarchical models ...
         }
      }

      // Parameter information vectors
      Rcpp::CharacterVector parNames;               // collected names of used hpar and par vectors (where size>0)
      Rcpp::LogicalVector parHyper;                 // if it is a hpar (TRUE) or par (FALSE) vector
      Rcpp::IntegerVector parSizes;                 // size (number of levels) in each hpar and par vector
      Rcpp::IntegerVector parEstFirst,parEstLast;   // First and last number in estimates vector (R coding from 1)
      Rcpp::IntegerVector parModelNr;               // which model-number the parameter vector comes from
      Rcpp::CharacterVector parModelFunc;           // The model-function name (fx, rn, rr, ...) the parameter vector comes from
      // idea to also add type, e.g. var, coeff, prob, to make easy selections in summary()
      Rcpp::IntegerVector parLogged;                // If the parameter will be logged
      Rcpp::CharacterVector parLoggedNames;         // Shortened list of names of logged parameters
      Rcpp::CharacterVector estimNames;             // Names of all 'estimates' (individual levels)
      
      collectParInfo(model, parNames, parHyper, parSizes, parEstFirst, parEstLast, parModelNr,
                     parModelFunc, parLogged, parLoggedNames, estimNames);
      size_t nEstimates = estimNames.size();

      // Check the chain settings and make list of output sample cycle-numbers.
      if (chain[0]==0 && chain[1]==0 && chain[2]==0) {  // chain was not set
         chain[0]=1100; chain[1]=100; chain[2]=10;
         Rcpp::Rcout << "Warning: chain was not set, running 1100 cycles but it may be too short for many analyses\n";
      }
      if (chain.size() != 3) throw (generalRbayzError("The chain settings do not have 3 elements"));
      if (chain[1] <= 0) throw (generalRbayzError("The chain length is zero or negative"));
      if (chain[2] < 0 || chain[3]<0 ) throw (generalRbayzError("The chain burnin or skip is negative"));
      if (chain[3] == 0) chain[3]=1;  // interpret skip zero as outputting all cycles
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
      Rcpp::rownames(loggedSamples) = Rcpp::as<Rcpp::CharacterVector>(outputCycleNumbers); 
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
               (Rcpp::Named("ModelNr")=parModelNr, Rcpp::Named("ModelTerm")=parModelFunc, 
                Rcpp::Named("Hyper")=parHyper,
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
      result.push_back(chain,"Chain");
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

// A function to insert one of the indepVarStr classes in the model vector.
// Note1: for MIXT it inserts two objects! But the actual variance object is the last one.
// Note2: if there is not varianceType set (0, no V= was found), the default IDEN is inserted;
// this part of code is only supposed to be called when a variance model is needed.
void insertIndepVar(std::vector<modelBase *> & model, dcModelTerm & modeldescr, modelBase * coeffmod) {
   if (modeldescr.varianceType==0 || modeldescr.varianceType==1) {
      // for now assume there are no multiple varianceNames and there are no parentheses with
      // parameters in the variance specifications
      if (modeldescr.varianceType==0 || modeldescr.varianceNames[0]=="IDEN") {
         model.push_back(new idenVarStr(modeldescr, coeffmod));
      }
      else if (modeldescr.varianceNames[0]=="MIXT") {
         // first insert mixture indicator class
         model.push_back(new mixtVarStr(modeldescr, coeffmod));
      }
      else {
         throw (generalRbayzError("Unknown variance term in model term "+modeldescr.funcName+": "+modeldescr.varianceNames[0]));
      }
   }
   else if (modeldescr.varianceType==2) {
      throw (generalRbayzError("Cannot use correlated variance structures for model term: "+modeldescr.funcName));
   }
   else {  // varianceType must be 3
      throw (generalRbayzError("Not yet ready to insert variance linear model for model term: "+modeldescr.funcName));
   }
}

void buildModel(std::vector<modelBase *> & model, dcModelTerm & modeldescr, modelBase * rmod)
{
   if (modeldescr.variableNames[0]=="1")
      model.push_back(new modelMean(modeldescr, rmod));
   else if (modeldescr.variableNames[0]=="0")
      return;
   else if (modeldescr.funcName=="fx" || modeldescr.funcName=="fixf")
      model.push_back(new modelFixf(modeldescr, rmod));
   else if (modeldescr.funcName=="rn" || modeldescr.funcName=="ranf") {
      if(modeldescr.varianceType==2) {
         // here need to add like for ranfi below
         model.push_back(new modelRanf_cor(modeldescr, rmod));
      }
      else {
         modelRanfi* tempptr = new modelRanfi(modeldescr, rmod);
         model.push_back(dynamic_cast<modelBase *>(tempptr));
         insertIndepVar(model, modeldescr, model.back());
         tempptr->varmodel = dynamic_cast<indepVarStr *>(model.back());
      }
   }
   else if (modeldescr.funcName=="rg" || modeldescr.funcName=="freg")
      model.push_back(new modelFreg(modeldescr, rmod));
   else if (modeldescr.funcName=="rr") {
      // For the moment only accept rr model with xxx/<matrix>
      if (modeldescr.variableNames.size()==2 && modeldescr.hierarchType==1 && modeldescr.variableTypes[1]==6) {
         modelRreg* tempptr = new modelRreg(modeldescr, rmod);
         model.push_back(dynamic_cast<modelBase *>(tempptr));
         insertIndepVar(model, modeldescr, model.back());
         tempptr->varmodel = dynamic_cast<indepVarStr *>(model.back());
      }
      else {
         throw (generalRbayzError("Cannot yet use rr() with something else than id/matrix"));
      }
      // Here should create an indepVarStr object with pointer in the Rreg (last model) object
   }
   else if (modeldescr.funcName=="smurf")
      throw (generalRbayzError("Aha! Caught you trying to smurf a variable into the model"));
   else
      throw (generalRbayzError("Unknown model (function) term: "+modeldescr.funcName));
}

/* old model build
void buildModelTerm(std::vector<modelBase *> & model, std::string modelTerm, Rcpp::DataFrame & data,
                    simpleMatrix & resid, size_t respnr) {
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
 */

void collectParInfo(std::vector<modelBase *> & model, Rcpp::CharacterVector & parNames,
                    Rcpp::LogicalVector & parHyper, Rcpp::IntegerVector & parSizes,
                    Rcpp::IntegerVector & parEstFirst, Rcpp::IntegerVector & parEstLast,
                    Rcpp::IntegerVector & parModelNr, Rcpp::CharacterVector & parModelFunc, Rcpp::IntegerVector & parLogged,
                    Rcpp::CharacterVector & parLoggedNames, Rcpp::CharacterVector & estimNames) {

   // First collect list of used hpar and par vectors, and set what model(-term) they belong to.
   for(size_t mt=0; mt<model.size(); mt++) {     // Done in two loops over model-terms, all hpar will come
      if ( model[mt]->hpar.nelem > 0 ) {         // first, then all par vectors
         parHyper.push_back(TRUE);
         parModelNr.push_back(mt);
         parModelFunc.push_back(model[mt]->fname);
         parSizes.push_back(model[mt]->hpar.nelem);
         parNames.push_back(model[mt]->hparName);
         if (model[mt]->logPars=="def") {           // default logging for hpar:
            if (model[mt]->hpar.nelem==1) {        // log hyper-parameter with 1 level
               parLogged.push_back(1);
               parLoggedNames.push_back(model[mt]->hparName);
            }
            else
               parLogged.push_back(0);
         }
         else if (model[mt]->logPars=="all") {
            parLogged.push_back(1);
            parLoggedNames.push_back(model[mt]->hparName);
         }
         else  // anything else than "def" and "all" is ignored and gives no logging of the parameter!
            parLogged.push_back(0);
      }
   }
   for(size_t mt=0; mt<model.size(); mt++) {
      if (model[mt]->par.nelem > 0 ) {
         parHyper.push_back(FALSE);
         parModelNr.push_back(mt);
         parModelFunc.push_back(model[mt]->fname);
         parSizes.push_back(model[mt]->par.nelem);
         parNames.push_back(model[mt]->parName);
         if (model[mt]->logPars=="def") {          // default logging for regular par:
            if (model[mt]->par.nelem==1) {        // log all parameters with 1 level. Could also consider to log
               parLogged.push_back(2);            // 2nd level of 2-level fixf parameters (code 3 ) ...
               parLoggedNames.push_back(model[mt]->parName);
            }
            else
               parLogged.push_back(0);
         }
         else if (model[mt]->logPars=="all") {
            parLogged.push_back(2);
            parLoggedNames.push_back(model[mt]->parName);
         }
         else  // anything else than "def" and "all" is ignored and gives no logging of the parameter!
            parLogged.push_back(0);
      }
   }
   
   // Fill names of all estimates (parameter%level), and fill parEstFirst, parEstLast which tell where each
   // parameter-vector can be found in the estimates list. This runs one loop, but that requires to set pointers to
   // par or hpar vectors (names, levels); alternative is to copy 2nd part of code two times in loops above....
   size_t nEstimates=1;  // to fill First and Last on index-base 1, after loop nEstimates
                         // will be 1 more than total levels
   {
      simpleDblVector *par_ptr;
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
         nEstimates += par_ptr->nelem;
         parEstLast.push_back(nEstimates-1);
         if (par_ptr->nelem==1) {  // for parameter-vectors size 1, estimate name is same as parameter name
            estimNames.push_back(*name_ptr);
         }
         else {                    // parameter-vector >1: the estimNames get appended %level
            if (par_ptr->nelem != label_ptr->size()) {
               // Here is a piece of code filling dummy labels (1,2,...) when the length of
               // the labels does not match the length of the parameter vector. This mostly happens
               // when label-vector is not filled, but also avoids conflict when the two don't match.
               for (size_t level=0; level < par_ptr->nelem; level++) {
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
      if(model[mt]->hpar.nelem > 0 ) {
         for(size_t i=0; i<model[mt]->hpar.nelem; i++) {
            postMean[k] += model[mt]->hpar[i];
            postSD[k] += (model[mt]->hpar[i])*(model[mt]->hpar[i]);
            k++;
         }
      }
   }
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->par.nelem > 0 ) {
         for(size_t i=0; i<model[mt]->par.nelem; i++) {
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
