#ifndef HTSVCF_H
#define HTSVCF_H

extern "C"{
    #include "query_regions.h"
}

class htsVCF {
public:
    htsVCF(std::string fname, std::vector<std::string> regs = {});
    ~htsVCF();
    char * fname() const;
    int nregs() const;
    char * line();
    const char **list_chroms();
    bool next();

private:
    htsFile *fp_;
    std::vector<char *> regions_;
    int nregs_;
    kstring_t str_;
    tbx_t *tbx_;
    int current_reg_;
    hts_itr_t *itr_ ;
};

#endif // HTSVCF_H