#include <string>
#include "stringStreamLite.h"

#ifndef _TOKENPOSITION_
#define _TOKENPOSITION_

template<typename scalar = int>
scalar tokenPosition(std::string s, std::string token) {
  stringStreamLite ss(s, ':');
  std::string tok;
  scalar k = 0;
  while( ss >> tok ) {
    if(tok == token) return k;
    k++;
  }
  return -1;
}

#endif
