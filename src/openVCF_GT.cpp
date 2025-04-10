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
Rcpp::XPtr<VCFReader<GT, int>> openVCF_GT(std::string filename, std::vector<std::string> regions) {
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    }    
    Rcpp::XPtr<VCFReader<GT, int>> pin(new VCFReader<GT, int>(filename, regions));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples_GT(Rcpp::XPtr<VCFReader<GT, int>> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
Rcpp::List getLine_GT(Rcpp::XPtr<VCFReader<GT, int>> pin) {

    Rcpp::IntegerVector G = Rcpp::wrap(pin->values);

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
bool next_GT(Rcpp::XPtr<VCFReader<GT, int>> pin) {
    return pin->next();
}
