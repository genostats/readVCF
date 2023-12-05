#ifndef HTSVCF_H
#define HTSVCF_H

extern "C"{
    #include "query_regions.h"
}

class htsVCF {
public:
    htsVCF(std::string fname, std::vector<std::string> regs = {});
    ~htsVCF();
    std::string fname() const;
    int nregs() const;
    char * line();
    std::vector<std::string> list_chroms();
    bool next();

private:
    htsFile *fp_;
    std::vector<std::string> regions_;
    int nregs_;
    kstring_t str_;
    tbx_t *tbx_;
    struct {
        int in_headers_;
        int current_reg_;
    } info_reg_;
    hts_itr_t *itr_ ;
};

#endif // HTSVCF_H