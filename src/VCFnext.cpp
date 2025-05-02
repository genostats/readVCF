#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include "stringStreamLite.h"
#include "VCFfield.h"

// [[Rcpp::export]]
bool VCFnext(Rcpp::List x) {
  Rcpp::XPtr<VCFReader> pin = x["xptr"];
  return pin->next();
}
