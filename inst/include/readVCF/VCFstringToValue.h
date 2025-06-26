#include "VCFfield.h"
#include "sto.h"
#include <cstdlib>
#include <cmath> // for nan

#ifndef _VCFstring2value_
#define _VCFstring2value_

template<VCFfield field, typename scalar> 
class VCFstringToValue;

/********************************* GT ******************************/

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

/********************************* DS ******************************/
// same but converts a DS string to a scalar (float or double)
template<typename scalar>
class VCFstringToValue<DS, scalar> {

public:

inline scalar operator()(const char * s, int le) {
  if(s[0] == '.') {
    if (sizeof(scalar) == sizeof(float)) return std::nanf("");
    if (sizeof(scalar) == sizeof(double)) return std::nan("");
    else throw std::runtime_error("Encountered a Nan in a type that doesn't support it");
  }
  // there will be another nan check in atof
  return std::atof(s);
}

inline scalar operator()(std::string str) {
  if (str == ".") {
    if (sizeof(scalar) == sizeof(float)) return std::nanf("");
    if (sizeof(scalar) == sizeof(double)) return std::nan("");
    else throw std::runtime_error("Encountered a Nan in a type that doesn't support it");
  }
  // there will be another nan check in sto
  return sto<scalar>(str);
}
};



#endif
