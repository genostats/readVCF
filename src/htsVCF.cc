#include <Rcpp.h>
#include "htsVCF.h"


extern "C"{
    #include "query_regions.h" //BEFORE htsFile used in line 6 
    htsFile *hts_opening(const char *, const char *);
}


// CONSTRUCTOR : opening the file and calculating nber of regions to go through 
htsVCF::htsVCF(std::string fname, std::vector<std::string> regs)
: info_reg_({1, 0}), nregs_(regs.size()), str_({0,0,0})
{
    fp_ = hts_opening(fname.c_str(),"r");
    if ( !fp_ ) 
    {
        Rcpp::stop("Couldn't open file\n");
    }

    enum htsExactFormat format = (&fp_->format)->format; // will always be vcf but to be sure
    if ((format != vcf))    
    {
        Rcpp::stop("Your file is not in a supported format, try a vcf file\n");
    }

    tbx_ = tbx_index_load3(fname.c_str(), NULL, HTS_IDX_SAVE_REMOTE);// before was a ternary 
    if (regs.empty()) {
        //getting an array of all chroms presents in the cvf file, and doing as if these were all the regions asked
        //so putting them in regions_ and updating nregs at the same time
        const char ** tmp = tbx_seqnames(tbx_, &nregs_);//gives back a const char**, but regions_ not const
        if (!tmp) throw std::runtime_error("Failed to look up regions in .tbi");
        regions_.reserve(nregs_);
        for (size_t i = 0; i < nregs_; ++i){   
            std::string reg(tmp[i]);
            regions_.push_back(reg);
        }
        free(tmp);
    } else {
        regions_ = regs;
    }
    itr_ = tbx_itr_querys(tbx_, regions_[info_reg_.current_reg_].c_str());// calls hts_itr_querys ret an itr or NULL if error
}



// DESTRUCTOR : deallocating some attributes, probably overkill in C++
htsVCF::~htsVCF()
{
   if (!regions_.empty()) { // happens if object never fully went through C°
        // TODO: rien ça marche
   }
    tbx_itr_destroy(itr_);
    tbx_destroy(tbx_);
    if (str_.s) free(str_.s);
    //now free fp : TODO : see if no double free !!!
    sam_hdr_destroy(fp_->bam_header);
    hts_idx_destroy(fp_->idx);
    hts_filter_free(fp_->filter);
    if (fp_->format.compression != no_compression) bgzf_close(fp_->fp.bgzf);
    free(fp_->fn);
    free(fp_->fn_aux);
    free(fp_->line.s);
    free(fp_);
}

// auto s = std::string(std::begin(txt_msg), std::find(std::begin(txt_msg), std::end(txt_msg), '\0')); 

// GETTERS
std::string htsVCF::fname() const { return std::string(fp_->fn, std::strlen(fp_->fn)); }
char * htsVCF::line() { return str_.s; }
int htsVCF::nregs() const { return nregs_; }



//Updating str_.s with the next line of data from the file.
bool htsVCF::next()
{
    if (info_reg_.current_reg_ >= nregs_) return false;
    int ret = 0;
    if (info_reg_.current_reg_ == 0 && info_reg_.in_headers_) // TODO : current_reg not really necessary in this check
    {
        if ((ret = hts_getline(fp_, '\n', &str_)) >= 0) // 2nd parameter unused, but must be '\n' (or KS_SEP_LINE) (else aborting)
        { // hts_getline updates fp_.lineno
            if ( !(!str_.l || str_.s[0]!= 35) ) {// 35 == conf.meta_char, == #
                return true;
            }
            else // if left the area with #
              info_reg_.in_headers_ = 0;
        }
        if (ret < -1) throw std::runtime_error("The reading of the file failed"); //-1 on end-of-file; <= -2 on error
    }
    while ( !itr_ ) { // could be false if region invalid
        if ( ++info_reg_.current_reg_ >= nregs_ )  {
            return false;// no more regions
        }
        else {
            Rcpp::Rcout << "Warning : " << regions_[info_reg_.current_reg_ ] << " seems to be an invalid region. \n";
            itr_ = tbx_itr_querys(tbx_, regions_[info_reg_.current_reg_].c_str());
        }
    }

    ret = tbx_itr_next(fp_, tbx_, itr_, &str_);// == hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx)) @return >= 0 on success, -1 when there is no more data, < -1 on error
    if (ret == -1) // when no more data, need to go to next reg
    {
        itr_ = NULL;
        return this->next();
    } 
    if (ret < -1)
    {
        throw std::runtime_error("An error occured while reading the line\n");
    }
    return true;

}

//Returns an array of all chromosomes in the target file.
std::vector<std::string> htsVCF::list_chroms()
{
    const char **seq = NULL;
    int i, nseq = 0;
    seq = tbx_seqnames(tbx_, &nseq); //very weird func who gives nbr of chroms to nseq and array of chroms to seq
    if (!seq) {
        free(seq);
        throw std::runtime_error("Failed to get list of sequence names");
    }
    std::vector<std::string> ret;
    for (int i = 0; i < nseq; ++i) {
        std::string str = seq[i];
        ret.push_back(str);// TODO : check if cast possible ? 
    }
    free(seq);
    return ret;

}
