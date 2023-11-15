#include <Rcpp.h>
#include "htsVCF.h"

// [[Rcpp::export]]
void test_pointer(std::string filename, std::vector<std::string> regions) {
  //Converting type std::string into char *;
  const char * fname = filename.c_str();
  std::vector<char*> regions_c; //normally compatible with char** in C 
  regions_c.reserve(regions.size());
  
  for(size_t i = 0; i < regions.size(); ++i)
    regions_c.push_back(const_cast<char*>(regions[i].c_str()));//un peu nasty...


  // creating htsVCF object
  htsVCF test_regs(filename, regions);

  // messing around with regions
  // (the vector could even be destroyed !)
  for(std::string & s : regions)
    s = "azerty";

  int i = 0;
  while (test_regs.next())
  {
    Rcpp::Rcout << test_regs.line() << '\n';
  }
}
