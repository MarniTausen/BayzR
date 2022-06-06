// correlVarStr.h --- Classes of correlated variance structures.
// Right now, there is only one type storing the eigendecomposition and in case of multiple kernels
// making the eigendecomposition of the kronecker product.
// Note: work in progress! Objects like ranf_cor do not yet interface with a generic correlVarStr
// like ranfi does, needs more work to develop.... this code still copy from indepVarStr....

#ifndef correlVarStr_h
#define correlVarStr_h

#include <Rcpp.h>
#include "modelVar.h"
#include "modelBase.h"
#include "dcModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"

class correlVarStr : public modelVar {
public:
   correlVarStr(dcModelTerm & modeldescr, modelBase* cm) : modelVar(modeldescr, cm) {
       coefmodel = cm;
       weights.initWith(cm->par.nelem, 1.0l);
   }
   simpleDblVector weights;
   modelBase* coefmodel;
};

class kernelVarStr : public correlVarStr {
public:
    kernelVarStr(dcModelTerm & modeldescr, modelBase* cm) : correlVarStr(modeldescr, cm) {
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

#endif /* correlVarStr_h */
