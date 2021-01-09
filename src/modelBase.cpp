#include "modelBase.h"

// modelBase constructor
modelBase::modelBase(dcModelTerm & modeldescr, modelBase * rmod)
                 : gprior(modeldescr.priorModel) {
   if (modeldescr.hierarchType == 0 || modeldescr.hierarchType == 1) {
      respModel = rmod;
      if (respModel != NULL) {
         resid = respModel->resid;
         residPrec = respModel->residPrec;
         Nresid = respModel->Nresid;
      }
   } else { // Hierarchical model:
            // don't know yet how to do here, there can be different cases,
            // but the rmod may not have a resid and residPrec vector
      respModel = rmod;
   }
   hierType = modeldescr.hierarchType;
   parName = modeldescr.allVariableNames;
   varNames = modeldescr.variableNames;
   varObjects = modeldescr.variableObjects;
   varType = modeldescr.variableTypes;
   parName = modeldescr.allVariableNames;
}
