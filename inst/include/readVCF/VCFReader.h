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
    std::vector<int> genos;
    VCFsnpInfo<int> snp;
    bool finished = false;
    VCFReader(std::string filename, std::vector<std::string> regions = {}) : in(filename,regions) {
        // skip VCF header
        while(in.next()) {
            char * li = in.line();
            if( (li[0] != '#') | (li[1] != '#') )
              break;
        }
        // on doit être sur la ligne qui contient les samples
        readVCFsamples(in.line(), samples);
        in.next();
        VCFlineGenotypes(in.line(), snp, genos);
    }
    bool next() {
        bool ret = in.next();
        genos.clear();
        if (ret) 
           VCFlineGenotypes(in.line(), snp, genos);
        else 
           finished = true;
        return ret;
    }
};

#endif
