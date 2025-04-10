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
Rcpp::XPtr<VCFReader<DS, double>> openVCF_DS(std::string filename, std::vector<std::string> regions) {
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    }    
    Rcpp::XPtr<VCFReader<DS, double>> pin(new VCFReader<DS, double>(filename, regions));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples_DS(Rcpp::XPtr<VCFReader<DS, double>> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
Rcpp::List getLine_DS(Rcpp::XPtr<VCFReader<DS, double>> pin) {

    Rcpp::NumericVector G = Rcpp::wrap(pin->values);

    return Rcpp::List::create(
        Rcpp::Named("chr") = pin->snp.chr,
        Rcpp::_("POS") = pin->snp.pos,
        Rcpp::_("ID") = pin->snp.id,
        Rcpp::_("REF") = pin->snp.ref,
        Rcpp::_("ALT") = pin->snp.alt,
        Rcpp::_("QUAL") = pin->snp.qual,
        Rcpp::_("FILTER") = pin->snp.filter,
        Rcpp::_("INFO") = pin->snp.info,
        Rcpp::_("genotypes") = G
    );
}

// [[Rcpp::export]]
bool next_DS(Rcpp::XPtr<VCFReader<DS, double>> pin) {
    return pin->next();
}
