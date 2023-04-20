#ifndef AHC_CPP_UTIL_H
#define AHC_CPP_UTIL_H

#include <stdexcept>
#include <string>

#define AHC_ASSERT_MSG(expression, message)\
do {\
    if (!(expression)) {\
        std::string msg;\
        msg += "Expression '" #expression "' is false. ";\
        msg += message;\
        throw std::runtime_error(msg);\
    }\
} while (false)

#define AHC_ASSERT(expression) AHC_ASSERT_MSG(expression, "Internal error.")

#endif // AHC_CPP_UTIL_H
