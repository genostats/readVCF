#include <Rcpp.h>
#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include "VCFfield.h"
#include "getField.h"


/**
 * @brief Write an integer dosage matrix (.dos16) from a VCF file.
 *
 * This function reads genotype dosage values from a VCF file, extracts the values
 * from a given FORMAT field (either "DS" or "GP"), and writes them into a binary
 * matrix file using 16-bit unsigned integers. It also generates a corresponding
 * `.bim` file with SNP metadata.
 *
 * Dosage values are scaled by 8192 and rounded to the nearest integer.
 * Missing values are encoded as 32768 (outside the 0–32767 range).
 *
 * @param filename Path to the input VCF file (bgzipped and tabix-indexed recommended).
 * @param newfile_name Output file prefix (without extension); will generate `.dos16` and `.bim`.
 * @param regions Vector of genomic regions to extract (e.g. {"19:320000-330000"}).
 *                If empty, the entire VCF is processed.
 * @param field The FORMAT field to extract dosage from. Supports "DS" or "GP". Defaults to "DS".
 *
 * @return An Rcpp List containing:
 *   - `samples`: a vector of sample IDs
 *   - `nbSNPs`: the number of variants successfully written
 *
 * @section Errors:
 * Throws std::runtime_error if:
 *   - the FORMAT field is not found in the VCF,
 *   - output files already exist (non-empty),
 *   - files cannot be created or written to.
 * 
 * @details
 * The `.dos16` file is written in row-major order (SNPs in rows, samples in columns),
 * as a sequence of `uint16_t` values. A matching `.bim` file is written with PLINK-style metadata.
 *
 * This function is intended to be called from R via Rcpp and does not overwrite existing files.
 */
// [[Rcpp::export]]
Rcpp::List writeDosageInteger(std::string filename, std::string newfile_name, std::vector<std::string> regions, std::string field = "DS") {
    // Initialize VCF reader with given regions
    VCFReader to_read(filename, regions);

    // Check that "DS" format is present in the VCF
    std::vector<std::string> formats = to_read.formats;
    if(std::find(formats.begin(), formats.end(), field) == formats.end())
        throw std::runtime_error("Field '" + field + "' not found in FORMATs of file " + filename);

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
            getField(to_read, field, val);// Read dosage values (DS field)
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
