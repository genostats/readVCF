#include <iostream>
#include <string>
#include <sstream>
#include "stringStreamLite.h"
#include "sto.h"
#ifndef TOKENATPOSITION
#define TOKENATPOSITION

template<typename T>
T tokenAtPosition(std::string & s, int pos) {
  stringStreamLite ss(s, ':');
  std::string token;
  for(int i = 0; i < pos && ss.next_token() != 0; i++) {}
  ss >> token;
  T r = sto<T>(token);
  return r;
}

#endif
