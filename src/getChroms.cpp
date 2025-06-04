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
Rcpp::CharacterVector getChroms(Rcpp::List x) {
  Rcpp::XPtr<VCFReader> pin = x["xptr"];
  return Rcpp::wrap(pin->in.list_chroms()); // ou sinon pin->in.list_chroms()
}