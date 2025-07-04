#include <iostream>
#include <cmath>
#include <cstring>        
#include <string_view>    

#include "../inst/include/readVCF/VCFstringToValue.h"


int main() {
  VCFstringToValue<GP, double> gp_parser;

  auto test = [&](const char* gp_str, double expected) {
    double result = gp_parser(gp_str, std::strlen(gp_str));
    if (std::isnan(expected)) {
      if (!std::isnan(result)) {
        std::cerr << "FAIL: " << gp_str << " → got " << result << ", expected NaN\n";
      } else {
        std::cout << "OK: " << gp_str << " → NaN\n";
      }
    } else {
      if (std::abs(result - expected) > 1e-6) {
        std::cerr << "FAIL: " << gp_str << " → got " << result << ", expected " << expected << "\n";
      } else {
        std::cout << "OK: " << gp_str << " → " << result << "\n";
      }
    }
  };

  test("0.0,0.0,1.0", 2.0);   // 0 + 2×1 = 2
  test("0.1,0.2,0.7", 1.6);   // 0.2 + 2×0.7 = 1.6
  test("1,0,0", 0.0);         // 0 + 2×0 = 0
  test("0.3,0.3,0.4", 1.1);   // 0.3 + 2×0.4 = 1.1
  test(".,.,.", std::nan("")); // NA

  return 0;
}
