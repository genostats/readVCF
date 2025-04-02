#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineGenotypes.h"


#ifndef __vcfreader__
#define __vcfreader__
class VCFReader {
    public:
    htsVCF in;
    std::vector<std::string> samples;
    char * line;
    VCFsnpInfo<int> snp;
    //std::string line;
    VCFReader(std::string filename, std::vector<std::string> regions) : in(filename, regions) {
        // skip VCF header
        while(in.next()) {
            char * li = in.line();
            if( (li[0] != '#') | (li[1] != '#') )
              break;
        }
        // on doit être sur la ligne qui contient les samples
        readVCFsamples(in.line(), samples);
        in.next();
        std::vector<int> genos;
        VCFlineGenotypes(in.line(), snp, genos);
        //SNPids.push_back(snp.id); // TODO : see if of use ?
        in.next();
    }

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
/*class VCFsnpInfo {
  public:
  int chr;
  int pos;
  std::string id;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
};
*/
};

#endif
