#ifndef VCFstring2geno
#define VCFstring2geno

// convert a VCF string to genotype 0, 1, 2, and 3 for NA
template<typename scalar>
inline scalar VCFstringToGeno(const char * s, int le) {
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

template<typename scalar>
inline scalar VCFstringToGeno(std::string str) {
  return VCFstringToGeno<scalar>(str.c_str(), str.length());
}

#endif
