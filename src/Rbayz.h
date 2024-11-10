//
//  Rbayz.h
//
// Some global variables used across constructors and functions:
//  - parList: is accessed in modelBase constructor (top of hierarchy) to collect vector of model parameters.
//  - Messages and needStop: can be used in any (helper) function finding errors. When functions not immediately
//    throw an exception, higher level code should check needStop and throw an exception.


#ifndef Rbayz_h
#define Rbayz_h

#include "parVector.h"

namespace Rbayz {
   extern std::vector<parVector**> parList;
   extern std::vector<std::string> Messages;
   extern bool needStop;
}

#endif /* Rbayz_h */
