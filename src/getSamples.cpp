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
Rcpp::CharacterVector getSamples(Rcpp::List x) {
  Rcpp::XPtr<VCFReader> pin = x["xptr"];
  return Rcpp::wrap( pin->samples ); 
}

