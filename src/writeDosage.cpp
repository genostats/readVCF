#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include "VCFfield.h"

// GOALS :
// fonction qui sauve une matrice de type "bigmemory"
// en données brutes => donc qui écrit 2 files,
// - un avec une extension .dosf (parce que des floats)
// - un avec une extension .dosf.bim, selon le format de PLINK
// => matrice créée à partir des dosages d'un fichier .vcf


// PARAMETERS :
// le path du file vcf
// le path du file à remplir (avec extension .dosf et .dosf.bim)
// POTENTIELLEMENT DES REGIONS !!! => si seulement les dosages 


// RETURN TYPE :
// tous les samples présents

// [[Rcpp::export]]
std::vector<std::string> writeDosage(std::string filename, std::string newfile_name, std::vector<std::string> regions) {
    VCFReader to_read(filename, regions);

    // check to_read.formats to make sure that there is dosage
    std::vector<std::string> formats =  to_read.formats;
    if (std::find(formats.begin(), formats.end(), "DS") == formats.end())
        throw std::runtime_error("No dosage format found for the file " + filename);

    /* Opening files for writing,
    my existence check will be simpler than snipsnop.
    If the file is empty, i will consider that I can overwrite it */

    // typeid(val[0]).name() => ret "f" for float. could be useful if templated
    // to have correct file extension maybe !
    // For now, hardcoded to float
    newfile_name += ".dosf";
    // app is for appending, binary is for writting raw data
    std::ofstream newfile(newfile_name, std::ios::binary | std::ios::app);
    if (!newfile)
        throw std::runtime_error("Error when attempting to create " + newfile_name);
    if (newfile.tellp()) // so if file size != 0 
        throw  std::runtime_error("Potentially overwriting an existing file \"" + newfile_name + "\" ... Aborting !");

    // opening .bim 
    std::string newfile_bim_name = newfile_name + ".bim";
    std::ofstream newfile_bim(newfile_bim_name, std::ios::binary | std::ios::app);
    if (!newfile_bim)
        throw std::runtime_error("Error when attempting to create " + newfile_bim_name);
    if (newfile_bim.tellp()) // so if file size != 0 
        throw  std::runtime_error("Potentially overwriting an existing file \"" + newfile_bim_name + "\" ... Aborting !");

    std::vector<float> val;

    while(to_read.next()) {
        bool ok = to_read.get<DS>(val);
        if (!ok)
            throw std::runtime_error("Something went wrong while reading the line");
        if (val.empty())
            throw std::runtime_error("Error fetching the dosage line");

        newfile.write(reinterpret_cast<const char *>(&val[0]), val.size() * sizeof(float));
        // Filling up the .bim with https://www.cog-genomics.org/plink/1.9/formats#bim
        newfile_bim << to_read.snpInfos.chr << "\t" << to_read.snpInfos.id << "\t0\t" // 0 hardcoded for pos in morgans
        << to_read.snpInfos.pos << "\t" << to_read.snpInfos.ref << "\t" << to_read.snpInfos.alt << "\n";
    }
    newfile.close();
    newfile_bim.close();

    return to_read.samples; 
}