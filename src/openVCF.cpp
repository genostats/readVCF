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
Rcpp::XPtr<VCFReader> openVCFregions(std::string filename, std::vector<std::string> regions) {
    /*
    for (auto reg : regions){ // will change regions to just an empty vector identifiable by .empty() in htsVCF c°
        if (reg.empty() && regions.size() == 1)  regions = {};
    } 
    */   
    Rcpp::XPtr<VCFReader> pin(new VCFReader(filename, regions));
    pin.attr("class") = "VCFReader";
    return pin;
}
