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
SEXP readVCFgenotypes(std::string filename, std::string regions = "") {
  // if (!(regions.empty()) || (filename.rfind(".gz") == filename.length() - 3))
  // {
  //   Rcpp::Rcout << "There are regions to map";
  //   //Converting type std::string into char *;
  //   const char * filename_c = filename.c_str();
  //   char ** arr_regions;
  //   printf("Splitting string %s into tokens:\n",regions);
  //   char *regions_c = &regions[0]; 
  //   arr_regions[0] = strtok(regions_c," ");
  //   int i = 1;
  //   while (arr_regions[i] != NULL)
  //   {
  //     arr_regions[i] = strtok (NULL, " ");
  //     printf ("this is first reg : %s\n",arr_regions[i]);
  //     i++;
  //   }
  //   interface(filename_c);
  // }
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
