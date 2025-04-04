#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include "VCFlineGenotypes.h" // used in next


// [[Rcpp::export]]
Rcpp::XPtr<VCFReader> openVCFregs(std::string filename, std::vector<std::string> regions) {
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    }    
    Rcpp::XPtr<VCFReader> pin(new VCFReader(filename, regions));
    return pin;
}

//[[Rcpp::export]]
Rcpp::XPtr<VCFReader> openVCFonly(std::string filename) {
    Rcpp::XPtr<VCFReader> pin(new VCFReader(filename));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples(Rcpp::XPtr<VCFReader> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
Rcpp::List getLine(Rcpp::XPtr<VCFReader> pin) {

    if (!(pin->in.line()) || pin->finished ) {return "No more lines to read";}
    
    Rcpp::IntegerVector G = Rcpp::wrap(pin->genos);

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
Rcpp::LogicalVector getNextLine(Rcpp::XPtr<VCFReader> pin) {
    Rcpp::LogicalVector ret = pin->in.next();
    if (ret[0]) VCFlineGenotypes(pin->in.line(), pin->snp, pin->genos);  
    else pin->finished++; 
    return ret;
}