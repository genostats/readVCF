#include <Rcpp.h>
#include "VCFReader.h"
#include "VCFlineValues.h"
#include "stringStreamLite.h"

using namespace Rcpp;

/**
 * @brief Read summary statistics from a VCF file and extract specified FORMAT fields.
 *
 * For each variant line in the VCF file, this function extracts standard VCF fields
 * and specific FORMAT fields (e.g., ES, SE, LP). It returns one column per sample and field,
 * with numeric values stored as float and text fields (e.g., ID) as string.
 *
 * If multiple samples are present, each sample gets its own set of columns (e.g., S1_AF, S2_AF).
 * Missing values (".") are correctly handled and converted to NA.
 *
 * @param vcf_file Path to the bgzipped VCF file.
 * @param fields Vector of FORMAT field names to extract (e.g., "ES", "SE", "LP").
 * @param regions Optional list of regions to restrict parsing (e.g., "1:1000-2000").
 *
 * @return A data.frame with standard VCF fields and per-sample FORMAT field columns.
 */
// [[Rcpp::export]]
DataFrame readSummaryStats(std::string vcf_file,
                           std::vector<std::string> fields,
                           Nullable<std::vector<std::string>> regions = R_NilValue) {

  // Convert regions (if provided) from Nullable to std::vector
  std::vector<std::string> regionsVec;
  if (regions.isNotNull()) {
    regionsVec = Rcpp::as<std::vector<std::string>>(regions);
  }

  // Open the VCF file using the VCFReader utility
  VCFReader reader(vcf_file, regionsVec);

  // Output vectors for fixed VCF columns
  std::vector<std::string> chrVec, idVec, refVec, altVec, qualVec, filterVec;
  std::vector<int> posVec;

  std::vector<std::string> sampleNames = reader.samples;
  size_t nSamples = sampleNames.size();

  // String-only fields (all others are parsed as float)
  std::set<std::string> stringFields = {"ID"};

  // Separate containers: float values and string values per sample × field
  std::vector<std::map<std::string, std::vector<float>>> sampleFieldValuesNum(nSamples);
  std::vector<std::map<std::string, std::vector<std::string>>> sampleFieldValuesStr(nSamples);

  // Initialize empty vectors for each requested field per sample
  for (size_t s = 0; s < nSamples; ++s) {
    for (const std::string &field : fields) {
      if (stringFields.count(field)) {
        sampleFieldValuesStr[s][field] = std::vector<std::string>();
      } else {
        sampleFieldValuesNum[s][field] = std::vector<float>();
      }
    }
  }

  // Process each line of the VCF
  while (reader.next()) {
    std::string line = reader.in.line();

    // Force VCFReader to parse current line
    std::vector<int> dummy;
    reader.get<VCFfield::GT, int>(dummy);
    const VCFsnpInfo<int> &snp = reader.snpInfos;

    // Split the raw line manually (more stable than relying on the reader)
    std::vector<std::string> tokens;
    size_t start = 0, end;
    while ((end = line.find('\t', start)) != std::string::npos) {
      tokens.push_back(line.substr(start, end - start));
      start = end + 1;
    }
    tokens.push_back(line.substr(start));

    if (tokens.size() < 9 + nSamples) {
      Rcpp::Rcout << "Invalid VCF line: " << line << std::endl;
      stop("VCF line has fewer than expected sample columns");
    }

    // Add standard VCF columns
    chrVec.push_back(tokens[0]);
    posVec.push_back(snp.pos);
    idVec.push_back(snp.id);
    refVec.push_back(snp.ref);
    altVec.push_back(snp.alt);
    qualVec.push_back(snp.qual);
    filterVec.push_back(snp.filter);

    // Parse FORMAT field definition (e.g., ES:SE:LP)
    std::vector<std::string> formatFields;
    {
      stringStreamLite parser(tokens[8], ':');
      while (parser) {
        std::string f;
        parser >> f;
        formatFields.push_back(f);
      }
    }

    // Process each sample column
    for (size_t s = 0; s < nSamples; ++s) {
      stringStreamLite valueParser(tokens[9 + s], ':');
      std::vector<std::string> values;
      while (valueParser) {
        std::string val;
        valueParser >> val;
        values.push_back(val);
      }

      for (size_t i = 0; i < formatFields.size(); ++i) {
        if (i < values.size()) {
          std::string fname = formatFields[i];
          if (stringFields.count(fname)) {
            sampleFieldValuesStr[s][fname].push_back(values[i]);
          } else {
            if (values[i] == "." || values[i].empty()) {
              sampleFieldValuesNum[s][fname].push_back(NA_REAL); // NA for float
            } else {
              sampleFieldValuesNum[s][fname].push_back(static_cast<float>(std::stof(values[i])));
            }
          }
        }
      }
    }
  }

  // Build final data.frame
  List df = List::create(
    _["chr"] = chrVec,
    _["pos"] = posVec,
    _["id"] = idVec,
    _["ref"] = refVec,
    _["alt"] = altVec,
    _["qual"] = qualVec,
    _["filter"] = filterVec
  );

  // Add one column per sample × field
  for (size_t s = 0; s < nSamples; ++s) {
    for (const std::string &field : fields) {
      std::string colName = sampleNames[s] + "_" + field;
      if (stringFields.count(field)) {
        df.push_back(sampleFieldValuesStr[s][field], colName);
      } else {
        df.push_back(NumericVector(sampleFieldValuesNum[s][field].begin(),
                                   sampleFieldValuesNum[s][field].end()),
                     colName);
      }
    }
  }

  df.attr("class") = "data.frame";
  df.attr("row.names") = IntegerVector::create(NA_INTEGER, -chrVec.size());
  return df;
}
