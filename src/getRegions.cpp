#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <stdint.h> //pour int64_t

// [[Rcpp::export]]
Rcpp::DataFrame getRegions_(Rcpp::List x) { // don't know which type would match best
  Rcpp::XPtr<VCFReader> pin = x["xptr"];
  
  interval_t *curr_interval = pin->in.list_intervals();

  if (!curr_interval) { // An exception should have been thrown in list_intervals !
    Rcpp::stop("Getting all regions failed");
  }
  std::vector<int> tid;
  std::vector<int> beg;  // should be a size_t but R won't convert those to int
  std::vector<int> end;

  while(curr_interval) {
    interval_t *to_free = curr_interval;
    tid.push_back(curr_interval->tid);
    beg.push_back(curr_interval->beg);
    end.push_back(curr_interval->end);

    curr_interval = curr_interval->next;
    free(to_free);
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create( Rcpp::Named("tid") = tid , Rcpp::_["beg"] = beg, Rcpp::_["end"] = end );
  return df;
}
