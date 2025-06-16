#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stringStreamLite.h"
#include "VCFfield.h"

// on renvoie un int, attention s'il y a plus de 2147483647 variants
// [[Rcpp::export]]
int countVariants_(std::string filename, std::vector<std::string> regions) {
  VCFReader reader(filename, regions);
  size_t n = 0;
  while(reader.next()) n++;
  if(n > 2147483647) Rcpp::stop("Integer too large for R data types");
  return (int) n;
}
