#ifndef GETFIELD_H
#define GETFIELD_H

#include <string>
#include <vector>
#include "VCFReader.h"
#include "VCFfield.h"

/**
 * @brief Extracts genotype information from a specified FORMAT field in a VCF line.
 *
 * This function dispatches to the appropriate template instantiation of `VCFReader::get<>`
 * based on the `field` string provided (e.g., "DS" or "GP"), and fills the `val` vector with the parsed values.
 *
 * @param reader Reference to an active VCFReader instance, positioned on a valid VCF line.
 * @param field  FORMAT field to extract ("DS" for dosage or "GP" for genotype probabilities).
 * @param val    Vector to store the resulting float values.
 *
 * @return True if the extraction succeeded; false otherwise.
 *
 * @section Errors:
 * Throws std::runtime_error If the field name is unsupported or invalid.
 */
inline bool getField(VCFReader & reader, const std::string & field, std::vector<float> & val) {
    if (field == "DS") return reader.get<DS>(val);
    else if (field == "GP") return reader.get<GP>(val);
    else throw std::runtime_error("Unsupported field: " + field);
}

#endif // GETFIELD_H
