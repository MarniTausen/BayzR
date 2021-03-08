#include "modelBase.h"

// modelBase constructor
modelBase::modelBase(dcModelTerm & modeldescr, modelBase * rmod)
{
   fname = modeldescr.funcName;
   hierType = modeldescr.hierarchType;
   parName = modeldescr.allVariableNames;
   varNames = modeldescr.variableNames;
   varObjects = modeldescr.variableObjects;
   varType = modeldescr.variableTypes;
   parName = modeldescr.allVariableNames;
}
