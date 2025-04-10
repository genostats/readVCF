#include "htsVCF.h"
#include "readVCFsamples.h"
#include "VCFsnpInfo.h"
#include "VCFlineValues.h"
#include "VCFfield.h"

#ifndef __vcfreader__
#define __vcfreader__
template<VCFfield field, typename scalar> 
class VCFReader {
    public:
    htsVCF in;
    std::vector<std::string> samples;
    std::vector<scalar> values;
    VCFsnpInfo<int> snp;
    VCFReader(std::string filename, std::vector<std::string> regions = {}) : in(filename,regions) {
        // skip VCF header
        while(in.next()) {
            char * li = in.line();
            if( (li[0] != '#') | (li[1] != '#') )
              break;
        }
        // on doit être sur la ligne qui contient les samples
        readVCFsamples(in.line(), samples);
        // attention on ne charge pas la première ligne maintenant. Il faut appeler next() avant de lire.
    }
    bool next() {
        bool ret = in.next();
        if(ret) {
          values.clear();
          VCFlineValues<field, scalar>(in.line(), snp, values);
        } 
        return ret;
    }
};

#endif
