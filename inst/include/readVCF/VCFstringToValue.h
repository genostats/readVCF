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
/********************************* GP ******************************/
// Converts a GP string like "0.1,0.2,0.7" to dosage: 0.2 + 2 * 0.7 = 1.6

template<typename scalar>
class VCFstringToValue<GP, scalar> {

public:

inline scalar operator()(const char * s, int le) {
  if (s[0] == '.') {
    if (sizeof(scalar) == sizeof(float)) return std::nanf("");
    if (sizeof(scalar) == sizeof(double)) return std::nan("");
    else throw std::runtime_error("Encountered a Nan in a type that doesn't support it");
  }

  // Parse manually (no need for stringStreamLite)
  scalar p0 = 0, p1 = 0, p2 = 0;
  const char* ptr = s;
  char* end;

  p0 = std::strtod(ptr, &end);
  if (*end != ',') goto parse_error;
  ptr = end + 1;

  p1 = std::strtod(ptr, &end);
  if (*end != ',') goto parse_error;
  ptr = end + 1;

  p2 = std::strtod(ptr, &end);
  // end should point to end of string or whitespace

  return p1 + 2.0 * p2;

parse_error:
  throw std::runtime_error(std::string("Invalid GP string: ") + s);
}

inline scalar operator()(std::string str) {
  return this->operator()(str.c_str(), str.length());
}
};



#endif
