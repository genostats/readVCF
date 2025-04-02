#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"



// [[Rcpp::export]]
Rcpp::XPtr<VCFReader> openVCF(std::string filename, std::vector<std::string> regions) {
    Rcpp::XPtr<VCFReader> pin(new VCFReader(filename, regions));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples(Rcpp::XPtr<VCFReader> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
Rcpp::String getLine(Rcpp::XPtr<VCFReader> pin) {
    return Rcpp::String(pin->line); // TODO : to fix print n'importe quoi
}

// [[Rcpp::export]]
Rcpp::LogicalVector getNextLine(Rcpp::XPtr<VCFReader> pin) {
    return Rcpp::wrap( pin->in.next() );
}
