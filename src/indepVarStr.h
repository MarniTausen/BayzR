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

    ~idenVarStr() {
        delete par;
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

/* diagVarStr is the model b~N(0,Ds^2), and D is diagonal, or on scalar level
   b_i ~ d_i s^2. It is the opposite of the weighted model where b_i ~ s^2 / w_i.
   diagVarStr is internally used for the rn_cor models running regression on eigenvectors,
   and D has the eigenvalues. 
*/
class diagVarStr : public indepVarStr {

public:

    diagVarStr(parsedModelTerm & modeldescr, parVector* coefpar, simpleDblVector & Ddiag)
            : indepVarStr(modeldescr, coefpar) {
        if(coefpar->nelem != Ddiag.nelem) {
            throw(generalRbayzError("ERROR dimension of DIAG does not fit random effect size"));
        }
        par = new parVector(modeldescr, 1.0l, "var");
        par->traced=1;
        par->varianceStruct="DIAG";
        // [ToDo]? this is copying the diagonal, but surely in the case when used for
        // eigenvector regression, this Ddiag are the eigenvalues that remain in the kernelMatrix
        // object, and it can be enough to 'wrap' the diag vector here around the other.
        // Probably, it can be assumed this is always OK when using this constructor where the
        // third arg is a simpleDblVector (then it must remain present somewhere else until cleanup).
        diag.initWith(Ddiag);
    }

    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefpar->nelem; k++)
         ssq += coefpar->val[k]*coefpar->val[k]/diag.data[k];
      par->val[0] = gprior.samplevar(ssq,coefpar->nelem);
      double invvar = 1.0l/par->val[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar / diag.data[k];
    }

    simpleDblVector diag;
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
