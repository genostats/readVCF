#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "VCFstringToValue.h"
#include "tokenPosition.h"
#include "tokenAtPosition.h"
#include "VCFlineValues.h"

// [[Rcpp::export]]
int test1(std::string s) {
  VCFstringToValue<GT, int> converter;
  return converter(s);
}

// [[Rcpp::export]]
int test2(std::string s, std::string tok) {
  return tokenPosition<int>(s, tok);
}

// [[Rcpp::export]]
std::string test3(std::string s, int pos) {
  return tokenAtPosition<std::string>(s, pos);
}


std::string testLine = "2	136401418	rs57232086	A	G	.	.	PR	GT	0/0	0/0	0/0	0/1	1/1";

// [[Rcpp::export]]
Rcpp::IntegerVector test4() {
  VCFsnpInfo<int> snp;
  std::vector<int> genos;
  VCFlineValues<GT>(testLine, snp, genos);
  return Rcpp::wrap(genos);
}

