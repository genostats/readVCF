#ifndef HTSVCF_H
#define HTSVCF_H
// TODO : add other functions called in c°
extern "C"{
    #include "query_regions.h"
}

class htsVCF {
public:
    htsVCF(const char *fname, char **regs = NULL, int nregs = 0);
    ~htsVCF();
    char * fname() const;
    int nregs() const;
    char * line();
    const char **list_chroms();
    bool next();

private:
    htsFile *fp_;
    char **regions_;
    int nregs_;
    kstring_t str_;
    tbx_t *tbx_;
    int current_reg_;
    hts_itr_t *itr_ ;
};

#endif // HTSVCF_H