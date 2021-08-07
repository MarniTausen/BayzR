//
//  rbayzExceptions.hpp
//  rbayz
//
//  Created by Luc Janss on 16/09/2018.
//  Copyright © 2018 Luc Janss. All rights reserved.
//

#ifndef rbayzExceptions_hpp
#define rbayzExceptions_hpp

#include <stdio.h>
#include <exception>
#include <stdexcept>
#include <string>

/* Now using one generalRbayzError, which is basically a runtime_error using the standard
 * runtime_error what() string to store and report an error message.
 * By making it a separate class, I can still separately catch errors from rbayz and other
 * c++ errors (bad alloc, underflow, overflow could be common c++ errors).
 */

class generalRbayzError : public std::runtime_error {
public:
   generalRbayzError(std::string s) : std::runtime_error(s) { }
};


#endif /* rbayzExceptions_hpp */
