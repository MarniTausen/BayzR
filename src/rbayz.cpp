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
void writeScreenLogHead(Rcpp::CharacterVector & parLoggedNames);
void writeScreenLog(size_t & cycle, Rcpp::NumericMatrix & loggedCumMeans, size_t save);
void collectPostStats(std::vector<modelBase *> & model, Rcpp::NumericVector & postMean,
                      Rcpp::NumericVector & postSD);
void collectResiduals(std::vector<modelBase *> &model, Rcpp::NumericMatrix &residuals);
void finishResiduals(std::vector<modelBase *> &model, Rcpp::NumericMatrix &residuals, size_t save);
void collectLoggedSamples(std::vector<modelBase *> & model, Rcpp::IntegerVector & parModelNr,
                          Rcpp::IntegerVector & parLogged, Rcpp::NumericMatrix & loggedSamples,
                          Rcpp::NumericMatrix & loggedCumMeans, size_t save);

// !Global variable! Vector to collect pointers to all par-vectors pointers in the model-objects.
// I thought this was too annoying to pass around between all constructors, now the modelBase class
// takes care that these pointers are collected.
std::vector<parVector**> parList;

// [[Rcpp::export]]
Rcpp::List rbayz_cpp(Rcpp::Formula modelFormula, SEXP VE, Rcpp::DataFrame inputData,
                     Rcpp::IntegerVector chain, int silent)
//                   note VE must be a string, it will be converted below
{

   // Vectors for messages defined outside the try-block, so it remains available in catch() part.
   Rcpp::CharacterVector Messages;
   std::string lastDone;

   try {     // normal execution builds a return list at the end of try{}, in case of
             // errors, catch() builds a return list with the messages vector defined above.

      // split the modelFormula in a list of response (LHS) and explanatory (RHS) terms.
      std::string formulaAsCppstring = convertFormula(modelFormula);
      removeSpaces(formulaAsCppstring);
      std::vector<std::string> modelTerms = splitModelTerms(formulaAsCppstring);
                                    // response will be modelTerms[0]
                                    // intercept is inserted as mn(0) or mn(1)
      lastDone="Parsing model";
      if (silent==9) Rcpp::Rcout << "Parsing model done\n";

      // build response object - including response variance structure
      std::string VEstr =  Rcpp::as<std::string>(VE);
      parsedModelTerm parsedResponseVariable(modelTerms[0], VEstr, inputData);
      // here still need to add selecting different response objects based on variance structure
      modelResp* modelR = new modelResp(parsedResponseVariable);   

      // Build vector of modelling objects from RHS terms (loop from term=1)
      std::vector<modelBase *> model;
      for(size_t term=1; term<modelTerms.size(); term++) {
         parsedModelTerm pmt(modelTerms[term], inputData);
         if(pmt.funcName=="mn") model.push_back(new modelMean(pmt, modelR));
         else if(pmt.funcName=="fx") model.push_back(new modelFixf(pmt, modelR));
         else if(pmt.funcName=="rn") {
            if(pmt.varianceStruct=="iden")
               model.push_back(new modelRanFacIden(pmt, modelR));
            else
               throw generalRbayzError("There is no class to model rn(...) with Variance structure " + pmt.varianceDescr);
            }
         else {
           throw generalRbayzError("Unknown model-function \'" + pmt.funcName + "\' at "+pmt.shortModelTerm);
         }
      }
      lastDone="Model building";
      if (silent==9) Rcpp::Rcout << "Model building done\n";

      // compute number of residuals (data points) and number of parameters in the model.
      // Response object is built first and parList[0] has residuals/fitted values.
      size_t nResiduals = *(parList[0])->nelem;
      size_t nParameters = 0;
      for(size_t i=1; i<parList.size(); i++)
         nParameters += *(parList[i])->nelem;
      if(!silent) Rcpp::Rcout << "Model built with " << nResiduals << " data points (incl NA) and "
                              << nParameters << " parameters\n";

      // Check the chain settings and find number of output samples by making list of output cycle-numbers.
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

      // Find the number of traced parameters and set-up matrix to store samples of traced parameters
      size_t nTracedParam=0;
      if(!silent) Rcpp::Rcout << "Saving full traces for:";
      for(size_t i=0; i<parList.size(); i++) {
         if( *(parList[i])->logged ) {
            nTracedParam += *(parList[i])->nelem;
            if(!silent) {
               for(size_t j=0; j< *(parList[i])->nelem; j++) Rcpp::Rcout << " " << *(parList[i])->parLabels[j];
            }
         }
      }
      if(!silent) Rcpp::Rcout << "\n";
      Rcpp::NumericMatrix tracedSamples(int(nSamples),int(nTracedParam));

      // to show convergence set "nShow" interval and make vector to hold previously shown solutions
      int nShow = chain[0]/10;
      Rcpp::NumericVector prevShowConv(int(nTracedParam), 1.0l);
      lastDone="Preparing to run MCMC";
      if (silent==9) Rcpp::Rcout << "Preparing to run MCMC done\n";

      // Run the MCMC chain
      // ------------------

      if(!silent) {
         Rcpp::Rcout << "Cycle, cumulative postMeans for traced parameters and [convergence]\n:";
      }
      for (size_t cycle=1, save=0, showconv=0; cycle <= chain[0]; cycle++) {
         for(size_t mt=0; mt<model.size(); mt++)
            model[mt]->sample();
         // At the 'skip' intervals (and after burn-in): 1) update posterior statistics using
         // collectStats(); 2) save MCMC samples for the 'traced' parameters 
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0) ) {  // save cycle
            for(size_t mt=0; mt<model.size(); mt++) model[mt]->prepForOutput();
            for(size_t i=0; i<parList.size(); i++) *(parList[i])->collectStats();
            for(size_t i=0, col=0; i<parList.size(); i++) {
               if( *(parList[i])->logged ) {
                  for(size_t j=0; j< *(parList[i])->nelem; j++) {
                     tracedSamples[save,col] = *(parList[i])->val[j]; col++;
               }
            }
            save++;  // save is counter for output (saved) cycles
         }
         if (cycle % nShow == 0 && !silent) {  // show convergence
            Rcpp::Rcout << cycle;
            double conv_change=0.0l, postmean;
            for(size_t i=0, col=0; i<parList.size(); i++) {
               if( *(parList[i])->logged ) {
                  for(size_t j=0; j< *(parList[i])->nelem; j++) {
                     postmean = *(parList[i])->postMean[j];
                     Rcpp::Rcout << " " << postmean;
                     // compute relative change compared to previously stored values, ...
                     conv_change = abs( (prevShowConv[col] - postmean) / prevShowConv[col] );
                     // and keep this postmean for the next convergence showing
                     prevShowConv[col] = postmean;
                     col++;
                  }
               }
            }
            conv_change /= double(nTracedParam);
            // do not show conv_change the first time
            if(showconv>0) Rcpp::Rcout << " [" << conv_change << "]";
            else showconv=1;
            Rcpp::Rcout << "\n";
         }
      }
      lastDone="Finished running MCMC";
      if (silent==9) Rcpp::Rcout << "Finished running MCMC\n";

      // Build tables to go in the output:
      // ---------------------------------

      // 1. "Parameter" information table (this is not including residuals)
      Rcpp::CharacterVector parNames, parModelFunc, parVariables;
      Rcpp::IntegerVector parSizes, parEstFirst, parEstLast, parlogged;
      for(size_t i=1; i<parList.size(); i++) {
         parNames.push_back(parList[i]->parName);
         parModelFunc.push_back(parList[i]->funcName);
         parVariables.push_back(parList[i]->)
      }


      Rcpp::DataFrame parInfo = Rcpp::DataFrame::create
               (Rcpp::Named("ModelNr")=parModelNr, Rcpp::Named("ModelTerm")=parModelFunc, 
                Rcpp::Named("Hyper")=parHyper,
                Rcpp::Named("Size")=parSizes, Rcpp::Named("EstStart")=parEstFirst,
                Rcpp::Named("EstEnd")=parEstLast, Rcpp::Named("Logged")=parLogged);
      parInfo.attr("row.names") = parNames;

      // 2. "Estimates" table with parName, parLabels, postMean and postSD
      Rcpp::CharacterVector allParNames(nParameters);
      Rcpp::CharacterVector allParLabels(nParameters);
      Rcpp::CharacterVector allPostMeans(nParameters);
      Rcpp::CharacterVector allPostSDs(nParameters);
      for(size_t i=1, row=0; i<parList.size(); i++) {
         for(size_t j=0; j< *(parList[i])->nelem; j++) {
            allParNames[row] = *(parList[i])->parName;
            allParLabels[row] = *(parList[i])->parLabel[j];
            allPostMeans[row] = *(parList[i])->postMean[j];
            allPostSDs[row] = sqrt( *(parList[i])->postVar[j] );
            row++;
         }
      }
      lastDone="Computing postMeans and PostSDs";
      Rcpp::DataFrame estimates = Rcpp::DataFrame::create
              (Rcpp::Named("Param")=allParNames, Rcpp::Named("Level")=allParamLabels,
               Rcpp::Named("postMean")=postMean, Rcpp::Named("postSD")=postSD);
//      estimates.attr("row.names") = estimNames;  // no longer rownames??
      lastDone="Setting up estimates dataframe";

      // 3. "Samples" table (this is the matrix tracedSamples with row and col-names added)
      Rcpp::CharacterVector sampleColNames;
      for(size_t i=0; i<parList.size(); i++) {
         if( *(parList[i])->logged ) {
            if(parList[i]->nelem==1) sampleColNames.push_back(parList[i]->parName);
            else {
               std::string s = parList[i]->parName;
               for(size_t j=0; j<parList[i]->nelem; j++)
                  sampleColNames.push_back(s + parList[i]->parLabels[j])
            }
         }
      }
      Rcpp::colnames(tracedSamples) = sampleColNames;
      Rcpp::rownames(tracedSamples) = Rcpp::as<Rcpp::CharacterVector>(outputCycleNumbers); 

      // 4. "Residuals" table: this is now postMean of fitted value (from par-vector in modelR) and residuals
      // computed from difference with data (Y, also still stored in modelR). If missing value, residual is NA.
      // Need to think about modifications for non-linear models, then residual may also need to be stored
      // and averaged because it is more difficult to computer from the Y.
      Rcpp::NumericVector fitval(nResiduals);
      Rcpp::NumericVector resid(nResiduals);
      for(size_t i=0; i<nResiduals; i++) {
         fitval[i] = modelR->par->postMean[i];
         if(modelR->missing[i]) resid[i] = NA_REAL;
         else resid[i] = modelR->Y->val[i] - modelR->par->postMean[i];
      }
      Rcpp::DataFrame residuals = Rcpp::DataFrame::create
              (Rcpp::Named("Fitval")=fitval, Rcpp::Named("Residual")=resid;     
      residuals.attr("row.names") = Rcpp::as<Rcpp::CharacterVector>(modelR->par->parLabels); 

      // Build the final return list
      Rcpp::List result = Rcpp::List::create();
      result.push_back(0,"nError");
      result.push_back(parInfo,"Parameters");
      result.push_back(loggedSamples,"Samples");
      result.push_back(estimates,"Estimates");
      result.push_back(residuals,"Residuals");
      result.push_back(chain,"Chain");
      lastDone="Filling return list";
      if (silent==9) Rcpp::Rcout << "Ready filling return list\n";

      // normal termination
      // ------------------
      return(result);

   } // end try{}

   // Catching (all kinds of) errors (and build return list with the error message)
   // ------------------------------

   catch (generalRbayzError &err) {
      errorMessages.push_back(err.what());
   }
   catch (std::exception &err) {
      std::string s = std::string(err.what()) + " after " + lastDone;
      errorMessages.push_back(s);
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
   if (modeldescr.varianceType==0) {   // this can happen, use rn() or rr() without V=
      model.push_back(new idenVarStr(modeldescr, coeffmod));
   }
   else if (modeldescr.varianceType==1) {
      // here insert linmod version
      throw (generalRbayzError("Not yet ready to insert variance linear model for model term: "+modeldescr.funcName));
   }
   else {  // type must be 2 and variance must be one of the INDEP structures
      if (modeldescr.varianceNames[0]=="IDEN") {
         model.push_back(new idenVarStr(modeldescr, coeffmod));
      }
      else if (modeldescr.varianceNames[0]=="MIXT") {
         // add here to insert mixture indicator class before the variance object
         model.push_back(new mixtVarStr(modeldescr, coeffmod));
      }
      else {
         throw (generalRbayzError("Unknown variance term in model term "+modeldescr.funcName+": "+modeldescr.varianceNames[0]));
      }
   }
}

// This is not yet used ...
//void insertCorrelVar(std::vector<modelBase *> & model, dcModelTerm & modeldescr, modelBase * coeffmod) {
//}

void buildModel(std::vector<modelBase *> & model, dcModelTerm & modeldescr, modelBase * rmod)
{
   if (modeldescr.variableNames[0]=="1")
      model.push_back(new modelMean(modeldescr, rmod));
   else if (modeldescr.variableNames[0]=="0")
      return;
   else if (modeldescr.funcName=="fx" || modeldescr.funcName=="fixf")
      model.push_back(new modelFixf(modeldescr, rmod));
   else if (modeldescr.funcName=="rn" || modeldescr.funcName=="ranf") {
      if(modeldescr.varianceStruct==2) {  // correlated var-structure
         // For the moment ranf_cor builds its own variance structure and only
         // accepts kernels, needs work to work the same as ranfi.
         model.push_back(new modelRanf_cor(modeldescr, rmod));
      }
      else {   // INDEP var-structure
         modelRanfi* tempptr = new modelRanfi(modeldescr, rmod);
         model.push_back(dynamic_cast<modelBase *>(tempptr));
         insertIndepVar(model, modeldescr, model.back());
         tempptr->varmodel = dynamic_cast<indepVarStr *>(model.back());
      }
   }
   else if (modeldescr.funcName=="rg" || modeldescr.funcName=="freg")
      model.push_back(new modelFreg(modeldescr, rmod));
   else if (modeldescr.funcName=="rr") {
      // For the moment only accept rr model with xxx/<matrix> and only INDEP var structures
      if(modeldescr.varianceStruct==2) {
         throw (generalRbayzError("Cannot use rr() with correlated variance structure"));
      }
      if (modeldescr.variableNames.size()==2 && modeldescr.hierarchType==0 && modeldescr.variableTypes[1]==6) {
         modelRreg* tempptr = new modelRreg(modeldescr, rmod);
         model.push_back(dynamic_cast<modelBase *>(tempptr));
         insertIndepVar(model, modeldescr, model.back());
         tempptr->varmodel = dynamic_cast<indepVarStr *>(model.back());
      }
      else {
         throw (generalRbayzError("Cannot yet use rr() with something else than id/matrix"));
      }
   }
   else if (modeldescr.funcName=="smurf")
      throw (generalRbayzError("Aha! Caught you trying to smurf a variable into the model"));
   else
      throw (generalRbayzError("Unknown model (function) term: "+modeldescr.funcName));
}

void collectParInfo(std::vector<modelBase *> & model, Rcpp::CharacterVector & parNames,
                    Rcpp::LogicalVector & parHyper, Rcpp::IntegerVector & parSizes,
                    Rcpp::IntegerVector & parEstFirst, Rcpp::IntegerVector & parEstLast,
                    Rcpp::IntegerVector & parModelNr, Rcpp::CharacterVector & parModelFunc, Rcpp::IntegerVector & parLogged,
                    Rcpp::CharacterVector & parLoggedNames, Rcpp::CharacterVector & estimNames) {

   // First collect list of used hpar and par vectors, and set what model(-term) they belong to.
   for(size_t mt=0; mt<model.size(); mt++) {     // Done in two loops over model-terms, all hpar will come
      if ( model[mt]->hpar.nelem > 0 && model[mt]->fname != "rp") {  // first, then all par vectors
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
      if (model[mt]->par.nelem > 0 && model[mt]->fname != "rp") {
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

void writeScreenLog(size_t & cycle, Rcpp::NumericMatrix & loggedCumMeans, size_t save) {
   // note: save-counter is already updated, so the last available save'd item is (save-1).
   // to compute change also save-2 is needed, and output can only be written if save >= 2.
   if (save < 2) return;
   double change=0.0l;
   Rcpp::Rcout << cycle;
   for(size_t i=0; i<loggedCumMeans.ncol(); i++) {
      change += abs( (loggedCumMeans(save-1,i) - loggedCumMeans(save-2,i) ) / loggedCumMeans(save-2,i) );
      Rcpp::Rcout << " " << loggedCumMeans(save-1,i) ;
   }
   Rcpp::Rcout << " [" << (100.0l*change/double(loggedCumMeans.ncol())) << "]\n";
}

void collectPostStats(std::vector<modelBase *> & model, Rcpp::NumericVector & postMean,
                      Rcpp::NumericVector & postSD) {
   size_t k=0;
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->hpar.nelem > 0 && model[mt]->fname != "rp") {
         for(size_t i=0; i<model[mt]->hpar.nelem; i++) {
            postMean[k] += model[mt]->hpar[i];
            postSD[k] += (model[mt]->hpar[i])*(model[mt]->hpar[i]);
            k++;
         }
      }
   }
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->par.nelem > 0 && model[mt]->fname != "rp") {
         for(size_t i=0; i<model[mt]->par.nelem; i++) {
            postMean[k] += model[mt]->par[i];
            postSD[k] += (model[mt]->par[i])*(model[mt]->par[i]);
            k++;
         }
      }
   }
}

// Collect sums of residuals and fitted values from the response model-objects
void collectResiduals(std::vector<modelBase *> &model, Rcpp::NumericMatrix &residuals) {
   size_t k=0;             // counter for response models
   size_t col1, col2;
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->fname=="rp") {
         col1 = k*2;       // column to store residuals for k'th response
         col2 = k*2 + 1;   // column to store fitted values for k'th response
         modelResp* tempptr = dynamic_cast<modelResp*>(model[mt]);
         for(size_t row=0; row<tempptr->par.nelem; row++) {
            residuals(row,col1) += tempptr->par.data[row];
         }
         for(size_t row=0; row<tempptr->par.nelem; row++) {
            residuals(row,col2) += tempptr->fit.data[row];
         }
         k++;
      }
   }
}

void finishResiduals(std::vector<modelBase *> &model, Rcpp::NumericMatrix &residuals,
                     size_t save) {
   size_t k=0;
   size_t col1, col2;
   double nsave = double(save);
   for(size_t mt=0; mt<model.size(); mt++) {
      if(model[mt]->fname=="rp") {
         col1 = k*2;
         col2 = k*2 + 1;
         modelResp* tempptr = dynamic_cast<modelResp*>(model[mt]);
         for(size_t row=0; row<tempptr->par.nelem; row++) {
            if (tempptr->missing[row])
               residuals(row,col1) = NA_REAL;
            else
               residuals(row,col1) /= nsave;
         }
         for(size_t row=0; row<tempptr->par.nelem; row++) {
            residuals(row,col2) /= nsave;
         }
         k++;
      }
   }
}

void collectLoggedSamples(std::vector<modelBase *> &model, Rcpp::IntegerVector & parModelNr,
                          Rcpp::IntegerVector & parLogged, Rcpp::NumericMatrix & loggedSamples,
                          Rcpp::NumericMatrix & loggedCumMeans, size_t save) {
   size_t k=0;
   double x=0.0l, update=0.0l;
   // save is counter for saved samples, it is counted in the loop in the main function and
   // should go from 0 to loggedSamples.nrow().
   if(save > (loggedSamples.nrow()-1)) {
      throw(generalRbayzError("Storage for logged samples out of bounds"));
   }
   for(size_t i=0; i<parLogged.size(); i++) {
      if (parLogged[i] > 0) {
         if (parLogged[i] == 1) x = model[parModelNr[i]]->hpar[0];
         else if (parLogged[i] == 2) x = model[parModelNr[i]]->par[0];
         else if (parLogged[i] == 3) x = model[parModelNr[i]]->par[1];
         loggedSamples(save,k) = x;
         if (save==0) {
            loggedCumMeans(save,k) = x;
         }
         else {
            update = (x-loggedCumMeans(save-1,k))/double(save+1);
            loggedCumMeans(save,k) = loggedCumMeans(save-1,k) + update;
         }
         k++;
      }
   }
}
