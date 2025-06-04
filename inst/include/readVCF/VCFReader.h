#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineValues.h"
#include "VCFfield.h"

#ifndef __vcfreader__
#define __vcfreader__
class VCFReader {
  bool ok = false; //end or error
  public:
  htsVCF in;
  std::vector<std::string> samples;
  // TODO : !! careful ! type of chromosomes are hardcoded
  VCFsnpInfo<int> snpInfos;
  std::vector<std::string> formats;
  VCFReader(std::string filename, std::vector<std::string> regions = {}) : in(filename,regions) {
    // go through VCF header
    while(in.next()) {
      char * li = in.line();
      if( (li[0] != '#') | (li[1] != '#') )
        break;
      // reading the formats.
      if(strncmp(li, "##FORMAT=<ID=", 13) == 0)
      {
        std::string id;
        std::istringstream li1(li);
        std::getline(li1, id, ',');
        formats.push_back(id.substr(13));
      }
    }
    // on doit être sur la ligne qui contient les samples
    readVCFsamples(in.line(), samples);
    // attention on ne charge pas la première ligne maintenant. Il faut appeler next() avant get().
  }

  bool next() {
    ok = in.next();
    return ok;
  }

  template<VCFfield field, typename scalar> 
  bool get(std::vector<scalar> & values) {
    if(ok) {
      values.clear();
      VCFlineValues<field, scalar>(in.line(), snpInfos, values);
      return true;
    } else {
      return false;
    }
  }

};

#endif
