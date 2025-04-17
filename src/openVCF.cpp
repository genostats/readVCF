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
Rcpp::XPtr<VCFReader> openVCF(std::string filename, std::vector<std::string> regions) {
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    }    
    Rcpp::XPtr<VCFReader> pin(new VCFReader(filename, regions));
    return pin;
}

// [[Rcpp::export]]
Rcpp::CharacterVector getSamples(Rcpp::XPtr<VCFReader> pin) {
   return Rcpp::wrap( pin->samples ); 
}

// [[Rcpp::export]]
bool VCFnext(Rcpp::XPtr<VCFReader> pin) {
    return pin->next();
}

// [[Rcpp::export]]
Rcpp::List getLine(Rcpp::XPtr<VCFReader> pin, std::string field) {

    Rcpp::List L;

    if(field == "GT") {
      std::vector<int> val;
      pin->get<GT>(val);
      L["values"] = Rcpp::wrap(val);
    } else if(field == "DS") {
      std::vector<double> val;
      pin->get<DS>(val);
      L["values"] = Rcpp::wrap(val);
    } else {
      Rcpp::warning("Unknown field");
    }

    L["CHR"] = pin->snpInfos.chr;
    L["POS"] = pin->snpInfos.pos;
    L["ID"] = pin->snpInfos.id;
    L["REF"] = pin->snpInfos.ref;
    L["ALT"] = pin->snpInfos.alt;
    L["QUAL"] = pin->snpInfos.qual;
    L["FILTER"] = pin->snpInfos.filter;
    L["INFO"] = pin->snpInfos.info;

    return L;
}


