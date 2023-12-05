#include <Rcpp.h>
#include "htsVCF.h"

// [[Rcpp::export]]
void test_pointer(std::string filename, std::vector<std::string> regions) {
  // creating htsVCF object
  htsVCF test_regs(filename, regions);

  // messing around with regions
  // (the vector could even be destroyed !)
  for(std::string & s : regions)
    s = "azerty";

  while (test_regs.next())
  {
    Rcpp::Rcout << test_regs.line() << '\n';
  }
}
