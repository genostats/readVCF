#include <Rcpp.h>
#include "htsVCF.h"

// [[Rcpp::export]]
int test_htsVCF( std::string filename, std::vector<std::string> regions) {
  
  //Converting type std::string into char *;
  {
     const char * fname = filename.c_str();
     std::vector<char*> regions_c; //normally compatible with char** in C 
     regions_c.reserve(regions.size());
   
     for(size_t i = 0; i < regions.size(); ++i)
       regions_c.push_back(const_cast<char*>(regions[i].c_str()));//un peu nasty...

    Rcpp::Rcout << " This is filename : " << fname << " this is the first region : " << regions_c[0] << "\n";
    htsVCF test_regs(fname, &regions_c[0], regions.size());
    Rcpp::Rcout << " This is nregs : " << test_regs.nregs()  << "\n";
    Rcpp::Rcout << " This is fname : " << test_regs.fname()  << "\n";
    int i = 0;
    const char **seq = test_regs.list_chroms();
    while (seq[i])
    {
      Rcpp::Rcout << "This is seq[" << i << "]  = " << seq[i] << "\n";
      i++;
    }
    while (test_regs.next())
    {
      Rcpp::Rcout << test_regs.line();
    }
  }
  return 0;
}
