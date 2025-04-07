#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineValues.h"
#include "htsVCF.h"


// version qui utilise hts lib / fichier indexé etc
extern "C" int interface(const char *filename);
//IMPORTANT DE GARDER L'INCLUDE APRES POUR QUE interface SOIT CONSIDEREE COMME DU C, SINON DEF de query_regions compilée comme du C++
#include "query_regions.h"

// [[Rcpp::export]]
SEXP readVCFgenotypes2(std::string filename, std::vector<std::string> regions) {
  
  htsVCF in(filename, regions);

  // first skip header
  while(in.next()) {
    char * li = in.line();
    if( (li[0] != '#') | (li[1] != '#') )
      break;
  }
  // on doit être sur la ligne qui contient les samples
  std::vector<std::string> samples;
  readVCFsamples(in.line(), samples);

  // maintenant on lit le reste du fichier
  VCFsnpInfo<int> snp;
  std::vector<int> genos;
  std::vector<std::string> SNPids;
  while(in.next()) {
    VCFlineValues<GT>(in.line(), snp, genos);
    SNPids.push_back(snp.id);
  }

  // Bricoler une matrice à partir d'un vecteur>
  Rcpp::IntegerVector G = Rcpp::wrap(genos);
  G.attr("dim") = Rcpp::Dimension( samples.size(), SNPids.size() );
  // lui ajouter dimnames
  Rcpp::List dimNames(2);
  dimNames[0] = Rcpp::wrap(samples);
  dimNames[1] = Rcpp::wrap(SNPids);
  G.attr("dimnames") = dimNames;

  return G;
  return Rcpp::wrap(samples);
}
