#include "VCFfield.h"

#ifndef _VCFstring2value_
#define _VCFstring2value_

template<VCFfield field, typename scalar> 
class VCFstringToValue;

// class with a () operator. An object of this class
// converts a GT string to genotype 0, 1, 2, and 3 for NA
template<typename scalar>
class VCFstringToValue<GT, scalar> {

public:

inline scalar operator()(const char * s, int le) {
  scalar g = 0;
  if(le == 3) { // cas diploide
    if(s[0] == '1') g++;
    if(s[2] == '1') g++;
    if(s[0] == '.' || s[2] == '.') g = 3; // missing value : NA
  } else if(le == 1) { // cas haploide
    if(s[0] == '1') g++;
    if(s[0] == '.') g = 3; // missing value : NA
  } else {
    g = 3;
  }
  return g;
}

inline scalar operator()(std::string str) {
  return this->operator()(str.c_str(), str.length());
}

};


#endif
