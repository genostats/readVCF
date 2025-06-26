#include <iostream>
#include <string>
#include <sstream>

#ifndef _STO_
#define _STO_

template<typename T>
T sto(const std::string & x);

template<>
inline double sto<double>(const std::string&x) {
  if (x == ".") return std::nan("");
  return std::stod(x);
}

template<>
inline float sto<float>(const std::string&x) {
  if (x == ".") return std::nanf("");
  return std::stof(x);
}

template<>
inline int sto<int>(const std::string&x) {
  return std::stoi(x);
}

template<>
inline std::string sto<std::string>(const std::string&x) {
  return x;
}

#endif

