#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineGenotypes.h"


extern "C" int interface(const char *filename);
//IMPORTANT DE GARDER L'INCLUDE APRES POUR QUE interface SOIT CONSIDEREE COMME DU C, SINON DEF de query_regions compilée comme du C++
#include "query_regions.h"

// [[Rcpp::export]]
SEXP readVCFgenotypes(std::string filename) {
  if (filename.rfind(".gz") == filename.length() - 3) //check if filename ends with .gz
  {
    Rcpp::Rcout << "Found a .gz file, will call my function later ! \n";
    //Converting type std::string into char *;
    const char * filename_c = filename.c_str();
    interface(filename_c);
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
