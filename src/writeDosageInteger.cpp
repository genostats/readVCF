#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include "VCFfield.h"

//' Write integer dosage matrix (.dos16) from VCF file
//'
//' This function reads a VCF file with DS (dosage) format and generates:
//' - A `.dos16` binary file containing scaled integer dosage values (uint16_t)
//' - A `.bim` file with SNP metadata (chr, id, pos, ref, alt)
//'
//' Dosage values are multiplied by 8192 and rounded to the nearest integer.
//' Missing values are encoded as 32768.
//'
//' @param filename The path to the input VCF file.
//' @param newfile_name The output file prefix (without extension).
//' @param regions A list of regions to read from the VCF.
//' @return A list with:
//'   - `samples`: vector of sample IDs
//'   - `nbSNPs`: number of SNPs written
//'
//' @throws std::runtime_error if input is invalid or output files exist already.
//'
//' @details The `.dos16` file is a binary matrix of dimension (number of SNPs × number of samples)
//'          stored in row-major order (each SNP as a row).
//'
//' @note This function appends no data. It will fail if the target output file already exists.
//'
// [[Rcpp::export]]
Rcpp::List writeDosageInteger(std::string filename, std::string newfile_name, std::vector<std::string> regions) {
    // Initialize VCF reader with given regions
    VCFReader to_read(filename, regions);

    // Check that "DS" format is present in the VCF
    std::vector<std::string> formats = to_read.formats;
    if(std::find(formats.begin(), formats.end(), "DS") == formats.end())
        throw std::runtime_error("No dosage format found for the file " + filename);

    // Open output files
    std::string out_dos16 = newfile_name + ".dos16";
    std::ofstream dos16(out_dos16, std::ios::binary | std::ios::app);
    if(!dos16)
        throw std::runtime_error("Error when creating " + out_dos16);
    if(dos16.tellp())
        throw std::runtime_error("Potential overwrite of \"" + out_dos16 + "\". Aborting.");

    std::string out_bim = newfile_name + ".bim";
    std::ofstream bim(out_bim, std::ios::binary | std::ios::app);
    if(!bim)
        throw std::runtime_error("Error when creating " + out_bim);
    if(bim.tellp())
        throw std::runtime_error("Potential overwrite of \"" + out_bim + "\". Aborting.");

    // Internal buffers
    std::vector<float> val;
    std::vector<uint16_t> val_int;
    int nbSNPs = 0;
    
    // Constants
    const float scale = 8192.0f;
    const uint16_t missing_code = 32768;

    // Main reading loop over SNPs
    while (to_read.next()) {
        const size_t n = to_read.samples.size();
        val.assign(n, std::nanf(""));  // Initialize all dosages to NaN

        try {
            to_read.get<DS>(val); // Read dosage values (DS field)
            if (to_read.snpInfos.id == "rs1578411") {
                Rcpp::Rcout << "RS1578411 values:\n";
                for (float v : val)
                    Rcpp::Rcout << v << " ";
                Rcpp::Rcout << "\n";
            }

        } catch (...) {
            // If DS not present for this variant, keep values as NaN
        }


        // Scale and convert to uint16_t
        val_int.clear();
        val_int.reserve(n);
        for (float d : val) {
            if (std::isnan(d)) {
                val_int.push_back(missing_code);
            } else {
                int scaled = static_cast<int>(std::round(d * scale));
                if (scaled < 0) scaled = 0;
                if (scaled > 32767) scaled = 32767;
                val_int.push_back(static_cast<uint16_t>(scaled));
            }
        }

        // Display NaN after filling val_int
        for (size_t i = 0; i < val.size(); ++i) {
            if (std::isnan(val[i])) {
                Rcpp::Rcout << "NaN detected at index " << i << " → written as " << val_int[i] << "\n";
            }
        }

        // Write to .dos16 (binary) and .bim (text metadata)
        dos16.write(reinterpret_cast<const char *>(&val_int[0]), val_int.size() * sizeof(uint16_t));
        bim << to_read.snpInfos.chr << "\t" << to_read.snpInfos.id << "\t0\t"
            << to_read.snpInfos.pos << "\t" << to_read.snpInfos.ref << "\t" << to_read.snpInfos.alt << "\n";
        nbSNPs++;
    }

    dos16.close();
    bim.close();

    // Return result to R
    return Rcpp::List::create(
        Rcpp::Named("samples") = to_read.samples,
        Rcpp::Named("nbSNPs") = nbSNPs
    );
}
