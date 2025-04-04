#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"

// [[Rcpp::export]]
int test_htsVCF( std::string filename, std::vector<std::string> regions) {

    htsVCF test_noregs(filename);
    Rcpp::Rcout << " This is nregs when no regs given : " << test_noregs.nregs()  << "\n";
    Rcpp::Rcout << " This is fname when no regs given : " << test_noregs.fname()  << "\n";

    int c = 0;
    while (test_noregs.next())
    {
      //Rcpp::Rcout << test_noregs.line() << '\n';
      //Rcpp::Rcout << "true" << '\n';
      c++;
    }
    Rcpp::Rcout << "This is the total of line went through when no regions were given : " << c << '\n';

    
  /*{
    htsVCF with_regs(filename, regions);
    Rcpp::Rcout << " This is nregs : " << with_regs.nregs()  << "\n";
    Rcpp::Rcout << " This is fname : " << with_regs.fname()  << "\n";
    int count_reg = 0;

    while (with_regs.next())
    {
      Rcpp::Rcout << with_regs.line() << '\n';
      if (with_regs.line()[0] != 35) count_reg++;
    }
    auto vec = with_regs.list_chroms();
    for (auto i : vec)
      Rcpp::Rcout << "i = " << i << std::endl;
    Rcpp::Rcout << "nbr of snips identified =" << count_reg << std::endl;
  }*/

  return 0;
}
