#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "stringStreamLite.h"
#include "tokenPosition.h"
#include "tokenAtPosition.h"
#include "VCFstringToGeno.h"
#include "VCFsnpInfo.h"

#ifndef _VCFlineGenotypes_
#define _VCFlineGenotypes_

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
// chrT = le type pour les chromosomes
// scalar = le type numeriques pour les genotypes
// !! la fonction push back les SNP dans le vecteur genotypes, sans se préoccuper des données
// !! qui peuvent déjà s'y trouver 
template<typename lineT, typename chrT, typename scalar>
void VCFlineGenotypes(lineT line, VCFsnpInfo<chrT> & snp, std::vector<scalar> & genotypes) {

  stringStreamLite li(line, 9); // 9 = tab separated
  std::string format;
  if(!(li >> snp.chr >> snp.pos >> snp.id >> snp.ref >> snp.alt >> snp.qual >> snp.filter >> snp.info >> format)) {
    throw std::runtime_error("VCF file format error");
  }
  
  int pos = tokenPosition(format, "GT");
  if(pos != -1) {
    std::string G;
    while(li >> G) {
      // conversion du token t1 en génotype
      std::string GT( tokenAtPosition<std::string>(G, pos) );
      scalar g = VCFstringToGeno<scalar>(GT);
      genotypes.push_back(g);
    }
  }
}


#endif
