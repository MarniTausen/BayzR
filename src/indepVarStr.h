// indepVarStr.h --- Classes of independent variance structures.
// These have a common "interface" of a vector of individual variances, which is
// defined in the base class and accessible through a pointer to base class.
// The actual form of this variance structure and implementation details for updating
// the variances can remain hidden for most (or all) uses of it.

#ifndef indepVarStr_h
#define indepVarStr_h

#include <Rcpp.h>
#include "modelVar.h"
#include "dcModelTerm.h"
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"

class indepVarStr : public modelVar {
public:
   indepVarStr(dcModelTerm & modeldescr) : modelVar(modeldescr) {
        // add constructor
   }
   simpleDblVector weights;
};

class idenVarStr : public indepVarStr {
public:
    idenVarStr(dcModelTerm & modeldescr) : indepVarStr((modeldescr) {
        // add constructor
    }
    void sample() {
        // add sample implementation
    }
};

class mixtVarStr : public indepVarStr {
public:
    mixtVarStr(dcModelTerm & modeldescr) : indepVarStr((modeldescr) {
        // add constructor
    }
    void sample() {
        // add sample implementation
    }
};


#endif /* indepVarStr_h */
