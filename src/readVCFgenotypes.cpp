#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineGenotypes.h"
#include "query_regions.h"


extern "C" {
  int interface_cpp(char *filename){
    return interface(filename);
  }
}

// [[Rcpp::export]]
SEXP readVCFgenotypes(std::string filename) {
  if (filename.rfind(".gz") == filename.length() - 3) //check if filename ends with .gz
  {
    Rcpp::Rcout << "Found a .gz file, will call my function later ! \n";
    //for the time being list_chroms = 0, so false
    interface(filename);
  }
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

  // Bricoler une matrice à partir d'un vecteur>
  Rcpp::IntegerVector G = Rcpp::wrap(genos);
  G.attr("dim") = Rcpp::Dimension( SNPids.size(), samples.size() );
  // lui ajouter dimnames
  Rcpp::List dimNames(2);
  dimNames[0] = Rcpp::wrap(SNPids);
  dimNames[1] = Rcpp::wrap(samples);
  G.attr("dimnames") = dimNames;

  return G;
}
