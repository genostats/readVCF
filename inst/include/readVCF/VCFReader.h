#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineValues.h"
#include "VCFfield.h"

#ifndef __vcfreader__
#define __vcfreader__
class VCFReader {
  bool ok = false; 
  public:
  htsVCF in;
  std::vector<std::string> samples;
  VCFsnpInfo<int> snpInfos;
  VCFReader(std::string filename, std::vector<std::string> regions = {}) : in(filename,regions) {
    // go through VCF header
    while(in.next()) {
      char * li = in.line();
      if( (li[0] != '#') | (li[1] != '#') )
        break;
      // si je suis ici c'est forcément que li commence par ##
      if( (li[2] == 'F') && (li[3] == 'O') && (li[4] == 'R') && (li[5] == 'M') && (li[6] == 'A') && (li[7] == 'T'))
      {
        std::string formatli(li);
        int id_pos = formatli.find("ID=");
        if ( id_pos != std::string::npos) snpInfos.format.push_back(formatli.substr(id_pos + 3, 2));
        //je suis partie du principe que format était forcément limité à 2 characters
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
