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
Rcpp::XPtr<VCFReader<GT, int>> openVCFregs(std::string filename, std::vector<std::string> regions) {
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    }    
    Rcpp::XPtr<VCFReader<GT, int>> pin(new VCFReader<GT, int>(filename, regions));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples(Rcpp::XPtr<VCFReader<GT, int>> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
Rcpp::List getLine(Rcpp::XPtr<VCFReader<GT, int>> pin) {

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
bool getNextLine(Rcpp::XPtr<VCFReader<GT, int>> pin) {
    return pin->next();
}
