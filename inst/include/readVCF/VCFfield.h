#ifndef _VCFfield_
#define _VCFfield_

#include <string>

enum VCFfield { GT, DS, GP };

constexpr std::string_view fieldName(VCFfield field) {
  switch(field) {
    case GT: return "GT";
    case DS: return "DS";
    case GP: return "GP";
    default: return "??";
  }
}
#endif
