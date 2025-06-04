#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
#include "VCFReader.h"


/* CAREFUL : this function will load the full vcf to memory 
because it will go through the whole file */

// [[Rcpp::export]]
SEXP readVCFgenotypes2(std::string filename, std::vector<std::string> regions) {
  VCFReader in(filename, regions);
  
  // maintenant on lit le reste du fichier
  std::vector<int> genos_full;
  std::vector<int> genos_line;
  std::vector<std::string> SNPids;
  while(in.next()) {
    if(in.get<GT>(genos_line)) {
      //adding the line (maybe there is a cleaner way to do this ?)
      //by appending the genos_full vec with the genos_line
      genos_full.insert(genos_full.end(), genos_line.begin(), genos_line.end());
    }
    else
      throw std::runtime_error("Error when retrieving the line !");
    SNPids.push_back(in.snpInfos.id);
  }


  if (!genos_full.size()) throw std::runtime_error("Problem with genos (is empty)");

  // Bricoler une matrice à partir d'un vecteur>
  Rcpp::IntegerVector G = Rcpp::wrap(genos_full);

  //std::cout << "G size = " << G.size() << "\n";

  std::vector<std::string> samples = in.samples;

  if (samples.empty() || SNPids.empty())  throw std::runtime_error("Problem with the dimensions");

  G.attr("dim") = Rcpp::Dimension( samples.size(), SNPids.size() );

  // lui ajouter dimnames
  Rcpp::List dimNames(2);
  dimNames[0] = Rcpp::wrap(samples);
  dimNames[1] = Rcpp::wrap(SNPids);
  G.attr("dimnames") = dimNames;

  return G;
}