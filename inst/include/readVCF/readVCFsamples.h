#include "stringStreamLite.h"
#include <stdexcept>

#ifndef _readVCFsamples_
#define _readVCFsamples_

/*
template<typename T1, typename T2>
void readVCFsamples(T1 line, T2 & samples);
*/

template<typename T1>
void readVCFsamples(T1 line, std::vector<std::string> & samples) {
  stringStreamLite li(line, 9); // 9 = tab separated
  std::string G;
  // 9 champs qui ne sont pas des samples
  for(int i = 0; i < 9; i++) {
    if(!(li >> G))
      throw std::runtime_error("VCF file format error");
  }
  // on push back les samples
  while(li >> G) {
    samples.push_back( G );
  }
}
#endif
