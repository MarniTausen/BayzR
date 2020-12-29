// indepVarStr.h --- Classes of independent variance structures.
// These have a common "interface" of a vector of individual variances, which is
// defined in the base class and accessible through a pointer to base class.
// The actual form of this variance structure and implementation details for updating
// the variances can remain hidden for most (or all) uses of it.

#ifndef indepVarStr_h
#define indepVarStr_h

#include <Rcpp.h>
#include "parseFunctions.h"
#include "rbayzExceptions.h"
#include "simpleVector.h"

class indepVarStr {
public:
   indepVarStr(std::string modelTerm);
   virtual ~indepVarStr() {}
   virtual void sample()=0;
   simpleDblVector weights, par;
   std::string parName;
   std::vector<std::string> parLabels;
};

class idenVarStr : public indepVarStr {
public:
   
}
class mixtVarStr : public indepVarStr {
public:
   
}

#endif /* indepVarStr_h */
