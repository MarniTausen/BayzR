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
#include "dcModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"

class indepVarStr : public modelVar {
public:
   indepVarStr(dcModelTerm & modeldescr, modelBase* cm) : modelVar(modeldescr, cm) {
       coefmodel = cm;
       weights.initWith(cm->par.nelem, 1.0l);
       fname = "vi";
   }
   simpleDblVector weights;
   modelBase* coefmodel;
};

class idenVarStr : public indepVarStr {
public:
    idenVarStr(dcModelTerm & modeldescr, modelBase* cm) : indepVarStr(modeldescr, cm) {
      par.initWith(1,1.0l);
      parName = "var." + parName;
    }
    void sample() {
      double ssq=0.0;
      for(size_t k=0; k < coefmodel->par.nelem; k++)
         ssq += coefmodel->par[k]*coefmodel->par[k];
      par[0] = gprior.samplevar(ssq,coefmodel->par.nelem);
      double invvar = 1.0l/par[0];
      for(size_t k=0; k < weights.nelem; k++) weights[k] = invvar;
    }
};

class mixtVarStr : public indepVarStr {
public:
    mixtVarStr(dcModelTerm & modeldescr, modelBase* cm) : indepVarStr(modeldescr, cm) {
        // add constructor
    }
    void sample() {
        // add sample implementation
    }
};

class loglinVarStr : public indepVarStr {
public:
    loglintVarStr(dcModelTerm & modeldescr, modelBase* cm) : indepVarStr(modeldescr, cm) {
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
