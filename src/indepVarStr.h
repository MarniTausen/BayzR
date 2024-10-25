// indepVarStr.h --- Classes of independent variance structures.
// These have a common "interface" of a vector of individual variances, which is
// defined in the base class and accessible through a pointer to base class.
// The actual form of this variance structure and implementation details for updating
// the variances can remain hidden for most (or all) uses of it.
// The log-linear variance model is also in this class, because it can interface in the
// same way with a vector of variances per random effect.

#ifndef indepVarStr_h
#define indepVarStr_h

#include <Rcpp.h>
#include "modelVar.h"
#include "modelBase.h"
#include "parsedModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"
#include "parVector.h"
#include <unistd.h>

class indepVarStr : public modelVar {
public:
    indepVarStr(parsedModelTerm & modeldescr, parVector* cpar) : modelVar(modeldescr), weights() {
       coefpar = cpar;
       weights.initWith(coefpar->nelem, 1.0l);
   }
   simpleDblVector weights;
};

class idenVarStr : public indepVarStr {

public:
    idenVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
       par = new parVector(modeldescr, 1.0l, "var");
       par->traced=1;
       par->varianceStruct="IDEN";
    }

    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefpar->nelem; k++)
         ssq += coefpar->val[k]*coefpar->val[k];
      par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
      double invvar = 1.0l/par->val[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
    }
};

// mixtVarStr is now standard 2-class mixture with pi0, pi1, v0, v1
class mixtVarStr : public indepVarStr {
public:
    mixtVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {

        // add constructor
    }

    void sample() {
        // add sample implementation
    }

//    getSetOptions(modeldescr.options["V"]); // not yet defined?
//    modelBVS *mixmod;
};

class loglinVarStr : public indepVarStr {
public:
    loglinVarStr(parsedModelTerm & modeldescr, parVector* coefpar) : indepVarStr(modeldescr, coefpar) {
        // add constructor
    }
    void sample() {
        // add sample implementation
    }
};

/* for weighted variance the SSQ statistic will be multiplied by weight
      for(size_t k=0; k< M->ncol; k++)
         ssq += par[k]*par[k]*weights[k];
*/

#endif /* indepVarStr_h */
