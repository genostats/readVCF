#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineGenotypes.h"

// [[Rcpp::export]]
SEXP readVCFgenotypes(std::string filename) {
  std::ifstream in(filename);   
  if(!in.good()) {
    Rcpp::stop("Couldn't open file\n");
  }
  std::string line;
  // first skip header
  while(std::getline(in, line)) {
    if(line.substr(0,2) != "##")
      break;
  }
  // on doit être sur la ligne qui contient les samples
  std::vector<std::string> samples;
  readVCFsamples(line, samples);

  // maintenant on lit le reste du fichier
  VCFsnpInfo<int> snp;
  std::vector<int> genos;
  std::vector<std::string> SNPids;
  while(std::getline(in, line)) {
    VCFlineGenotypes(line, snp, genos);
    SNPids.push_back(snp.id);
  }

  // Bricoler une matrice à partir d'un vecteur
  Rcpp::IntegerVector G = Rcpp::wrap(genos);
  G.attr("dim") = Rcpp::Dimension( SNPids.size(), samples.size() );
  // lui ajouter dimnames
  Rcpp::List dimNames(2);
  dimNames[0] = Rcpp::wrap(SNPids);
  dimNames[1] = Rcpp::wrap(samples);
  G.attr("dimnames") = dimNames;

  return G;
}
