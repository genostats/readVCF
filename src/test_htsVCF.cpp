#include <Rcpp.h>
#include "htsVCF.h"

// [[Rcpp::export]]
int test_htsVCF( std::string filename, std::vector<std::string> regions) {

    htsVCF test_noregs(filename);
    Rcpp::Rcout << " This is nregs when no regs given : " << test_noregs.nregs()  << "\n";
    Rcpp::Rcout << " This is fname when no regs given : " << test_noregs.fname()  << "\n";

    while (test_noregs.next())
    {
      Rcpp::Rcout << test_noregs.line() << '\n';
    }

  {
    htsVCF with_regs(filename, regions);
    Rcpp::Rcout << " This is nregs : " << with_regs.nregs()  << "\n";
    Rcpp::Rcout << " This is fname : " << with_regs.fname()  << "\n";
    while (with_regs.next())
    {
      Rcpp::Rcout << with_regs.line() << '\n';
    }

  }

  return 0;
}
