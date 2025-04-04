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
Rcpp::String getLine(Rcpp::XPtr<VCFReader> pin) {
    std::string output = std::to_string(pin->snp.chr);
    std::string s{'\t'};
    std::string nl{'\n'};
    std::string pos =  std::to_string(pin->snp.pos);
    output = output + s + pos + s + pin->snp.id + s + pin->snp.ref + s + pin->snp.alt + s + pin->snp.qual + s + pin->snp.filter + s + pin->snp.info + nl;

    if (pin->in.line()) {
        std::string line(pin->in.line());
        //add line if safe to do so
        output = output + line + nl;
    }
    return Rcpp::String(output);
}

// [[Rcpp::export]]
Rcpp::LogicalVector getNextLine(Rcpp::XPtr<VCFReader> pin) {
    Rcpp::LogicalVector ret = pin->in.next();
    if (ret) VCFlineGenotypes(pin->in.line(), pin->snp, pin->genos);        
    return ret;
}
