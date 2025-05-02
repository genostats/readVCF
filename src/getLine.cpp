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
Rcpp::List getLine(Rcpp::List x, std::string field) {
    Rcpp::XPtr<VCFReader> pin = x["xptr"];
    
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
      L["values"] = R_NilValue;
      Rcpp::warning("Unknown format field");
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


