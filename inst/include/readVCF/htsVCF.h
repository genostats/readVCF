#include <string>
#include <vector>
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
    interval_t *list_chroms();
    bool next();

private:
    struct {
        int in_headers_;
        int current_reg_;
    } info_reg_;
    int nregs_;
    kstring_t str_;
    htsFile *fp_;
    std::vector<std::string> regions_;
    tbx_t *tbx_;
    hts_itr_t *itr_ ;
    bool idx_file_;
};
#endif // HTSVCF_H
