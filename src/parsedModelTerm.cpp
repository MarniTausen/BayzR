//  parsedModelTerm.cpp

#include <Rcpp.h>
#include <vector>
#include <string>
#include "parsedModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"

// parseModelTerm_step1: splits a model-term in 3 strings according to possible syntaxes
//  (1)   funcname(variableString,optionString)
//  (2)   funcname(variableString)
//  (3)   variableString
// Return is vector of 3 strings for funcname, variableString, optionString, some can
// be empty ("") if not available from syntaxes (2) and (3).
std::vector<std::string> parseModelTerm_step1(std::string mt) {

   std::vector<std::string> result(3);
   size_t pos1, pos2;
   bool has_funcname=false;

   // determine if there are parenthesis - so there should be a funcname
   if ( (pos1 = mt.find('(')) != std::string::npos)
      has_funcname=true;

   if (has_funcname && mt[mt.size()-1] != ')') {
      throw(generalRbayzError("No closing parenthesis in model-term: "+mt));
   }

   // get funcName and reset pos1 to where variable-list should start
   if ( has_funcname ) {
      result[0]=mt.substr(0, pos1);
      pos1++;
   }
   else {
      result[0]="";
      pos1=0;
   }

   // Determine if and where is ')' or ',' - this marks end of variableString with first 2
   // syntax patterns.
   // Checks: with a funcname, there must be ')' or ',' (but is already checked that there is ')')
   //         without a funcname, there cannot be ')' or ','
   pos2 = mt.find_first_of("),",pos1);
   if (!has_funcname && pos2 != std::string::npos) {     // wrong syntax
      if(mt[pos2]==')')
         throw(generalRbayzError("Unexpected closing parenthesis in response or model-term: "+mt));
      else
         throw(generalRbayzError("Unexpected comma in response or model-term: "+mt));
   }

   // 2. variableString
   size_t retrieve_length;
   if(pos2==std::string::npos)
      retrieve_length=pos2;            // 3rd syntax pattern: get substr to end
   else
      retrieve_length=pos2-pos1;       // 1st or 2nd syntax pattern: get substr to comma or )
   result[1]=mt.substr(pos1,retrieve_length);

   // 3. optionString: if pos2 is on a comma there are options upto closing parenthesis
   if(mt[pos2]==',') {
      retrieve_length = mt.size()-pos2-2;
      result[2]=mt.substr(pos2,retrieve_length);
   }
   else
      result[2]="";

   return result;

}

// parseModelTerm_step2: splitting / interpreting variables, options, etc.
// This one is defined as a member function to fill object member variables
parsedModelTerm::parseModelTerm_step2(std::string funcName, std::string variableString, std:string optionString) {

   if(variableString.length()<=12) {
     if(optionString=="") shortModelTerm=mt;
     else shortModelTerm=funcName+"("+variableString+",...)";
   }
   else {
     shortModelTerm=funcName+"("+variableString.substr(0,12)+"...)";
   }

   // Analyse and split variableString
   pos1 = variableString.find(':');
   pos2 = variableString.find('/');
   pos3 = variableString.find('|');
   if(pos2==std::string::npos && pos3==std::string::npos)            // A, A:B, A:B:C
      variablePattern="factors";
   else if (pos3!=std::string::npos && pos2==std::string::npos       // A|B, A|B:C
                  && (pos1==std::string::npos || pos2>pos3))         // but not A:B|C
      variablePattern="nestreg";
   else if (pos2!=std::string::npos && pos1==std::string::npos       // A/B but no other
                  && pos3==std::string::npos)                        // patterns allowed with /
      variablePattern="rrcovars";
   else {
      std::string s="Cannot interpret/use variable specification \'"+variableString+
                    "\' in model-term: "+mt;
      throw(generalRbayzError(s));
   }
   variableNames = splitString(variableString,":|/");

   // For every variable get an RObject pointing to it (whether it is from the data frame
   // or from R environment), and also store the type (factor, numeric, etc.) of the variable.
   for(size_t i=0; i<variableNames.size(); i++) {
      if (variableNames[i]=="1" || variableNames[i]=="0") {
         variableObjects.push_back(R_NilValue);
         variableTypes.push_back(0);
      }
      else {
         variableObjects.push_back(getVariableObject(d,variableNames[i]));
             // getVariableObject searches both the data frame 'd' and the R environment
         if(variableObjects.back() != R_NilValue)
            variableTypes.push_back(getVariableType(variableObjects.back()));
         else {
            throw generalRbayzError("Variable not found in data frame or R environment: "+variableNames[i]);
         }
      }
   }

   // split and store the options.
   // The options cannot simply be separated on commas, because there can be commas inside options
   // like: V=XX[a=1,b=2],W=YY[...]
   // Approach is therefore to scan character by character, check open & close brackets and compute
   // 'open_close_brack_balance'. A proper splitting comma is a comma where open_close_brack_balance==0. 
   // The split options are directly checked for being one of the know options (so far only V= and prior=)
   // and store the string after '=' in varianceDescr and priormodDescr.
   if (optionString!="") {
      pos1=-1;                          // start of first option, it will move up 1 at the start of the loop
      pos2=0;                           // will move to comma after first option, or end of string
      pos3=optionString.size()-1;       // position of last character in optionString
      int open_close_brack_balance=0;
      std::string tmpstring;
      do {
         pos1++;                        // for proper continuation: after processing an option,
                                        // pos1 will be left standing on the splitting comma.
         while( !(open_close_brack_balance==0 && optionString[pos2]==',') && pos2<pos3) {
            if(optionString[pos2]=='(' || optionString[pos2]=='[')
               open_close_brack_balance++;
            if(optionString[pos2]==')' || optionString[pos2]==']')
               open_close_brack_balance--;
            pos2++;
         }
         if(pos2==pos3)         // pos2 on the last character
            tmpstring=optionString.substr(pos1,(pos2-pos1+1));
         else                   // pos2 is after the last character (of piece to extract)
            tmpstring=optionString.substr(pos1,(pos2-pos1));
         if(tmpstring.substr(0,2)=="V=") {
            varianceDescr=tmpstring.substr(2,(tmpstring.size()-2));
         }
         else if (tmpstring.substr(0,6)=="prior=") {
            priormodDescr=tmpstring.substr(6,(tmpstring.size()-6));
         }
         else {
               std::string s="Unknown option in " + funcName + "(" + variableString +
                "...) : " + tmpstring;
                throw(generalRbayzError(s));
         }
         if(pos2==',') pos1=pos2;  // the while will continue for a next option
      }
      while (optionString[pos1]==',');
   }

   // Split and analyse the variance description. This writes in varianceStruct a string
   // that allows to select the right object class in main.
   if (varianceDescr=="") {
      varianceStruct="iden";
   }
   else {
      if (varianceDescr[0]=='~') {
         varianceStruct="llin";
      }
      else {    // all other cases should be structure keywords and kernels separated by stars
         std::vector<std::string> varianceElements = splitString(varianceDescr,"*");
         // Here there could be a way to allow fixing variances by detecting if the last
         // element is a numerical value, then store and remove that last element.
         // Split every variance element in a name and a parameter-part
         for(size_t i=0; i<varianceElements.size(); i++) {
            size_t bracket = varianceElements[i].find_first_of("([");
            size_t len_tot = varianceElements[i].length();
            std::string name;
            std::string params;
            if (bracket == std::string::npos) {
               name = varianceElements[i];
               params = "";
            }
            else {
               size_t closeBrack = findClosingBrack(varianceElements[i], bracket);
               if(closeBrack != (len_tot-1) ) {
                  throw generalRbayzError("Unbalanced parentheses in: "+varianceElements[i]);
               }
               name = varianceElements[i].substr(0,bracket);
               params = varianceElements[i].substr(bracket+1,(len_tot-bracket-2));
            }
            varianceNames.push_back(name);
            varianceParams.push_back(params);
         }
         // fill the varianceObjects vector (pointers to R objects)
         // all the elements with known bayz variance structure keywords get R_NilValue
         for(size_t i=0; i<varianceElements.size(); i++) {
            std::string name=varianceNames[i];
            if(name=="IDEN" || name=="MIXT" || name=="VCOV" || name=="WGHT" || name=="DIAG"
                    || name=="LLIN") {
               varianceObjects.push_back(R_NilValue);
            }
            else {
               varianceObjects.push_back(getVariableObject(d,name));
               if(varianceObjects.back() == R_NilValue) {
                  throw generalRbayzError("Variance/kernel object not found in the R Environment: "
                       +varianceNames[i]);
               }
            }
         }
         // and finally analyse the names and set varianceStruct to indicate what variance
         // class to use for the modelling ... but there are multiple terms and combinations ...
         if(varianceElements.size() == 1) {
            std::string name1=varianceNames[0];
            if(name1=="IDEN") varianceStruct="iden";
            else if(name1=="MIXT") varianceStruct="mixt";
            else if(name1=="WGHT") varianceStruct="wght";
            else if(name1=="DIAG") varianceStruct="diag";
            else if(name1=="VCOV") varianceStruct="vcov"; // only 1 VCOV unusual?
            else if(name1=="LLIN") varianceStruct="llin"; // needs fixing to handle it the same as V=~...
            else
               varianceStruct="kernel";
         }
         else {
            throw generalRbayzError("Cannot handle this variance pattern: "+varianceDescr);
         }
      }
   }

}

// constructor for handling response term with separate variance description
parsedModelTerm::parsedModelTerm(std::string mt, std::string VEdescr, Rcpp::DataFrame &d);
{
   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);
   // for now not accepting functions on response, but it could be extended here to
   // allow e.g. log(Y), probit(Y), etc
   if(parse_step1[0]!="") throw(generalRbayzError("Unexpected function on response term: "+mt));
   // the next one could just be a message, but not easy here to get things in the messages list
   if(parse_step1[2]!="") throw(generalRbayzError("Unexpected options retrieved from response term: "+mt));

   parseModelTerm_step2("", parse_step1[1], VEdescr);
   
}

// constructor for handling RHS model terms
parsedModelTerm::parsedModelTerm(std::string mt, Rcpp::DataFrame &d)
{
   size_t pos1, pos2, pos3;

   std::vector<std::string> parse_step1 = parseModelTerm_step1(mt);

   parseModelTerm_step2(parse_step1[0], parse_step1[1], parse_step1[2]);

}

