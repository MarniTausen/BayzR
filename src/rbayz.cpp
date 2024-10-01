#include <vector>
#include <string>
#include <Rcpp.h>
#include <cmath>
#include "parseFunctions.h"
#include "parsedModelTerm.h"
#include "modelBase.h"
#include "modelResp.h"
#include "modelMean.h"
#include "modelFixf.h"
#include "modelRanfi.h"
#include "modelFreg.h"
#include "modelRreg.h"
#include "modelRanfc.h"
#include "rbayzExceptions.h"
#include "simpleMatrix.h"
#include "simpleVector.h"
#include "modelVar.h"
#include "indepVarStr.h"
#include <unistd.h>

// [[Rcpp::plugins("cpp11")]]

// A few global variables are used! Otherwise it would need to pass these around in many functions
// and constructors.
// parList: is accessed in modelBase constructor (top of hierarchy) to collect vector of model parameters
// Messages and needStop: can be used in any (helper) function finding errors. When functions not immediately
// throw an exception, higher level code should check needStop and throw an exception.
// Note: the global variables appear to persist between R calls, they need to be cleared at the start of main().
std::vector<parVector**> parList;
std::vector<std::string> Messages;
bool needStop=false;

// [[Rcpp::export]]
Rcpp::List rbayz_cpp(Rcpp::Formula modelFormula, SEXP VE, Rcpp::DataFrame inputData,
                     Rcpp::IntegerVector chain, int verbose)
//                   note VE must be a string, it will be converted below
{

   // clearing global variables
   parList.clear();
   Messages.clear();
   needStop=false;

   // rbayz retains a small string of last executed code that is sometimes added in errors
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
      if (verbose > 1) Rcpp::Rcout << "Parsing model done\n";

      // build response object - including response variance structure
      std::string VEstr =  Rcpp::as<std::string>(VE);
      parsedModelTerm parsedResponseVariable(modelTerms[0], VEstr, inputData);
      // here still need to add selecting different response objects based on variance structure
      modelResp* modelR = new modelResp(parsedResponseVariable);   
      if (verbose > 1) Rcpp::Rcout << "Response model-object done\n";

      // Build vector of modelling objects from RHS terms (loop from term=1)
      std::vector<modelBase *> model;
      for(size_t term=1; term<modelTerms.size(); term++) {
         parsedModelTerm pmt(modelTerms[term], inputData);
         if(pmt.funcName=="mn") model.push_back(new modelMean(pmt, modelR));
         else if(pmt.funcName=="fx") model.push_back(new modelFixf(pmt, modelR));
         else if(pmt.funcName=="rn") {
            if(pmt.varianceStruct=="iden" || pmt.varianceStruct=="notgiven")
               model.push_back(new modelRanFacIden(pmt, modelR));
            else
               throw generalRbayzError("There is no class to model rn(...) with Variance structure " + pmt.options["V"]);
         }
         else if (pmt.funcName=="rr") {
            if(pmt.varianceStruct=="iden" || pmt.varianceStruct=="notgiven")
               model.push_back(new modelRregIden(pmt, modelR));
            else
               throw generalRbayzError("There is no class to model rr(...) with Variance structure " + pmt.options["V"]);
         }
         else {
           throw generalRbayzError("Unknown model-function \'" + pmt.funcName + "\' at "+pmt.shortModelTerm);
         }
      } // end for(term ...) to build model
      lastDone="Model building";
      if (verbose>1) Rcpp::Rcout << "Model building done\n";

      // compute number of residuals (data points) and number of parameters in the model.
      // Response object is built first and parList[0] has residuals/fitted values.
      size_t nResiduals = (*(parList[0]))->nelem;
      size_t nParameters = 0;
      for(size_t i=1; i<parList.size(); i++)
         nParameters += (*(parList[i]))->nelem;
      if(verbose > 0) Rcpp::Rcout << "Model built with " << nResiduals << " data points (incl NA) and "
                              << nParameters << " parameters\n";

      if(verbose > 2) {
         Rcpp::Rcout << "Model-object overview (#, Name, Size, Logged, first Labels) after model building:\n";
         for(size_t i=0; i<parList.size(); i++) {
            Rcpp::Rcout << i << " " << (*(parList[i]))->Name << " " << (*(parList[i]))->nelem << 
                    " " << (*(parList[i]))->logged;
            Rcpp::Rcout << " " << (*(parList[i]))->Labels[0];
            if((*(parList[i]))->nelem > 1) Rcpp::Rcout << " " << (*(parList[i]))->Labels[1];
            if((*(parList[i]))->nelem > 2) Rcpp::Rcout << " ...";
            Rcpp::Rcout << "\n";
         }
      }

      // Make parameter name disambiguation - this skips residuals (parList[0])
      {
         Rcpp::CharacterVector parNames;
         for(size_t i=1; i<parList.size(); i++)
            parNames.push_back((*parList[i])->Name);
         Rcpp::IntegerVector name_matches = Rcpp::match(parNames, parNames);
         std::vector<size_t> numb_matches(name_matches.size(),0);
         bool found_duplicates = FALSE;
         for(int i=0; i < name_matches.size(); i++) {
            if(name_matches[i] != (i+1))  {  // the i'th name matches (name_matches[i]-1)'th name in the list
               if(numb_matches[name_matches[i]-1] == 0) {
                  parNames[name_matches[i]-1] += ".1";
                  numb_matches[name_matches[i]-1]++; 
               }
               numb_matches[name_matches[i]-1]++;
               parNames[i] += "." + std::to_string(numb_matches[name_matches[i]-1]);
            }
            found_duplicates=TRUE;
         }
         if(found_duplicates) {  // copy new names back in parList->Names
            for(int i=0; i<parNames.size(); i++)
               (*(parList[i+1]))->Name = parNames[i];
         }
      }

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
      for (int cycle=1; cycle <= chain[0]; cycle++) {         // exactly the same loop and condition to make
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0))  // output samples below
            outputCycleNumbers.push_back(cycle);
      }
      size_t nSamples = outputCycleNumbers.size();
      if (nSamples==0) throw (generalRbayzError("The chain settings do not make any output"));

      // Find the number of traced parameters and set-up matrix to store samples of traced parameters
      size_t nTracedParam=0;
      if(verbose>0) Rcpp::Rcout << "Saving full traces for:";
      for(size_t i=0; i<parList.size(); i++) {
         if( (*(parList[i]))->logged ) {
            nTracedParam += (*(parList[i]))->nelem;
            if(verbose>0) {
               for(size_t j=0; j< (*(parList[i]))->nelem; j++) Rcpp::Rcout << " " << (*(parList[i]))->Labels[j];
            }
         }
      }
      if(verbose>0) Rcpp::Rcout << "\n";
      Rcpp::NumericMatrix tracedSamples(nSamples,nTracedParam);

      // to show convergence set "nShow" interval and make vector to hold previously shown solutions
      int nShow = chain[0]/10;
      if( nShow < 1) nShow=1;
      Rcpp::NumericVector prevShowConv(int(nTracedParam), 1.0l);
      lastDone="Preparing to run MCMC";
      if (verbose>1) Rcpp::Rcout << "Preparing to run MCMC done\n";

      // Run the MCMC chain
      // ------------------

      if(verbose>0) {
         Rcpp::Rcout << "Cycle, cumulative postMeans for traced parameters and [convergence]\n";
      }
      for (int cycle=1, save=0, showconv=0; cycle <= chain[0]; cycle++) {
         modelR->sample();
         for(size_t mt=0; mt<model.size(); mt++) model[mt]->sample();
/*         for(size_t i=0; i<parList.size(); i++)
            Rcpp::Rcout << " " << (*(parList[i]))->val[0];
         Rcpp::Rcout << "\n";*/
         // At the 'skip' intervals (and after burn-in): 1) update posterior statistics using
         // collectStats(); 2) save MCMC samples for the 'traced' parameters 
         if ( (cycle > chain[1]) && (cycle % chain[2] == 0) ) {  // save cycle
            for(size_t mt=0; mt<model.size(); mt++) model[mt]->prepForOutput();
            for(size_t i=0; i<parList.size(); i++) (*(parList[i]))->collectStats();
            for(size_t i=0, col=0; i<parList.size(); i++) {
               if( (*(parList[i]))->logged ) {
                  for(size_t j=0; j< (*(parList[i]))->nelem; j++) {
                     tracedSamples(save,col) = (*(parList[i]))->val[j]; col++;
                  }
               }
            }
            save++;  // save is counter for output (saved) cycles
         }
         if (cycle % nShow == 0 && verbose>0) {  // show convergence
            Rcpp::Rcout << cycle;
            double conv_change=0.0l, postmean;
            for(size_t i=0, col=0; i<parList.size(); i++) {
               if( (*(parList[i]))->logged ) {
                  for(size_t j=0; j< (*(parList[i]))->nelem; j++) {
                     postmean = (*(parList[i]))->postMean[j];
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
      } // end for(cycle ...)
      lastDone="Finished running MCMC";
      if (verbose>1) Rcpp::Rcout << "Finished running MCMC\n";

      // Build tables to go in the output:
      // ---------------------------------

      // 1. "Parameter" information table
      Rcpp::CharacterVector parNames, parModelFunc, parVariables, parVarStruct;
      Rcpp::IntegerVector parSizes, parEstFirst, parEstLast, parLogged;
      for(size_t i=0, row=1; i<parList.size(); i++) {
         parNames.push_back((*(parList[i]))->Name);
         parModelFunc.push_back((*(parList[i]))->modelFunction);
         parVariables.push_back((*(parList[i]))->variables);
         parVarStruct.push_back((*(parList[i]))->varianceStruct);
         parSizes.push_back((*(parList[i]))->nelem);
         parEstFirst.push_back(row);
         parEstLast.push_back(row+(*(parList[i]))->nelem-1);
         row += (*(parList[i]))->nelem;
         parLogged.push_back((*(parList[i]))->logged);
      }
      Rcpp::DataFrame parInfo = Rcpp::DataFrame::create
               (Rcpp::Named("ModelTerm")=parModelFunc, Rcpp::Named("Variables")=parVariables, 
                Rcpp::Named("Param")=parNames,Rcpp::Named("Variance")=parVarStruct,
                Rcpp::Named("Size")=parSizes, Rcpp::Named("Start")=parEstFirst,
                Rcpp::Named("End")=parEstLast, Rcpp::Named("Logged")=parLogged);
//      parInfo.attr("row.names") = parNames;

      // 2. "Estimates" table with parName, parLabels, postMean and postSD
      Rcpp::CharacterVector allParNames(nParameters);
      Rcpp::CharacterVector allParLabels(nParameters);
      Rcpp::CharacterVector allPostMeans(nParameters);
      Rcpp::CharacterVector allPostSDs(nParameters);
      for(size_t i=1, row=0; i < parList.size(); i++) {
         for(size_t j=0; j< (*(parList[i]))->nelem; j++) {
            allParNames[row] = (*(parList[i]))->Name;
            allParLabels[row] = (*(parList[i]))->Labels[j];
            allPostMeans[row] = (*(parList[i]))->postMean[j];
            allPostSDs[row] = sqrt( (*(parList[i]))->postVar[j] );
            row++;
         }
      }
      lastDone="Computing postMeans and PostSDs";
      Rcpp::DataFrame estimates = Rcpp::DataFrame::create
              (Rcpp::Named("Param")=allParNames, Rcpp::Named("Level")=allParLabels,
               Rcpp::Named("postMean")=allPostMeans, Rcpp::Named("postSD")=allPostSDs);
//      estimates.attr("row.names") = estimNames;  // no longer rownames??
      lastDone="Setting up estimates dataframe";

      // 3. "Samples" table (this is the matrix tracedSamples with row and col-names added)
      Rcpp::CharacterVector sampleRowNames = Rcpp::wrap(modelR->par->Labels);
      Rcpp::CharacterVector sampleColNames;
      for(size_t i=0; i<parList.size(); i++) {
         if( (*(parList[i]))->logged ) {
            if( (*(parList[i]))->nelem==1) sampleColNames.push_back( (*(parList[i]))->Name);
            else {
               std::string s = (*(parList[i]))->Name;
               for(size_t j=0; j< (*(parList[i]))->nelem; j++)
                  sampleColNames.push_back(s + (*(parList[i]))->Labels[j]);
            }
         }
      }
      Rcpp::colnames(tracedSamples) = sampleColNames;
      Rcpp::rownames(tracedSamples) = Rcpp::as<Rcpp::CharacterVector>(outputCycleNumbers); 
      /* I couldn't get this colnames() and rownames() working, it gives a compiler error that the Rcpp::NumericMatrix
         can't be conveted to SEXP object - but online examples show this should work ...
         */
      // work-around to set row and colnames of a matrix: add "dimnames" attribute as a list of 2
//      Rcpp::List tracedSamplesNames = Rcpp::List::create(sampleRowNames,sampleColNames);
//      tracedSamples.attr("dimnames") = tracedSamplesNames;

      // 4. "Residuals" table: this is now postMean of fitted value (from par-vector in modelR) and residuals
      // computed from difference with data (Y, also still stored in modelR). If missing value, residual is NA.
      // Need to think about modifications for non-linear models, then residual may also need to be stored
      // and averaged because it is more difficult to computer from the Y.
      Rcpp::NumericVector fitval(nResiduals);
      Rcpp::NumericVector resid(nResiduals);
      for(size_t i=0; i<nResiduals; i++) {
         fitval[i] = modelR->par->postMean[i];
         if(modelR->missing[i]) resid[i] = NA_REAL;
         else resid[i] = modelR->Y.data[i] - modelR->par->postMean[i];
      }
      Rcpp::DataFrame residuals = Rcpp::DataFrame::create
              (Rcpp::Named("Fitval")=fitval, Rcpp::Named("Residual")=resid);
      Rcpp::CharacterVector residRowNames = Rcpp::wrap(modelR->par->Labels);
      residuals.attr("row.names") = residRowNames;

      // Build the final return list
      Rcpp::List result = Rcpp::List::create();
      result.push_back(0,"nError");
      result.push_back(parInfo,"Parameters");
      result.push_back(tracedSamples,"Samples");
      result.push_back(estimates,"Estimates");
      result.push_back(residuals,"Residuals");
      result.push_back(chain,"Chain");
      lastDone="Filling return list";
      if (verbose>1) Rcpp::Rcout << "Ready filling return list\n";

      // normal termination
      // ------------------
      return(result);

   } // end try{}

   // Catching (all kinds of) errors (and build return list with the error message)
   // ------------------------------

   catch (generalRbayzError &err) {
      Messages.push_back(err.what());
   }
   catch (std::exception &err) {
      std::string s = std::string(err.what()) + " after " + lastDone;
      Messages.push_back(s);
   }
   catch (...) {
      Messages.push_back("An unknown error occured in bayz after: "+lastDone);
   }

   // Note: program flow only comes here in case of errors (normal return is above);
   // build a return list that only has the error messages list.
   Rcpp::List result = Rcpp::List::create();
   result.push_back(Messages.size(),"nError");
   result.push_back(Messages,"Errors");
   Rcpp::Rcout << "Bayz finished with error: " << Messages[Messages.size()-1] << std::endl;
   Rcpp::Rcout << "There may be more messages or errors - use summary() or check <output>$Errors to see all" << std::endl;
   return(result);

}  // end rbayz_cpp main function

