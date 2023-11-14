#include <Rcpp.h>
#include "htsVCF.h"

// TODO: test with 



extern "C"{
    #include "query_regions.h" //BEFORE htsFile used in line 6 
    htsFile *hts_opening(const char *, const char *);
}


// CONSTRUCTOR : opening the file and calculating nber of regions to go through 
htsVCF::htsVCF(const char *fname,  char ** regs, int nregs)
: current_reg_(0), nregs_(nregs), str_({0,0,0}), regions_(nullptr)
{
    fp_ = hts_opening(fname,"r");
    if ( !fp_ ) 
    {
        Rcpp::stop("Couldn't open file\n");
    }

    enum htsExactFormat format = (&fp_->format)->format; // will always be vcf but to be sure
    if ((format != vcf)) 
    {
        Rcpp::stop("Your file is not in a supported format, try a vcf file\n");
    }

    tbx_ = tbx_index_load3(fname, NULL, HTS_IDX_SAVE_REMOTE);// before was a ternary 
    if (!regs || !(*regs)) {
        Rcpp::Rcout << "No region specified so will be going through the entire file\n If you wish filter it, please enter a region like so : Chrom_id:First_Index-Last_Index \n";
        //getting an array of all chroms presents in the cvf file, and doing as if these were all the regions asked
        //so putting them in regions_ and updating nregs at the same time
        regions_ = const_cast<char **>(tbx_seqnames(tbx_, &nregs_));//gives back a const char**, but regions_ not const
        if (!regions_)
        throw std::runtime_error("Failed to look up regions in .tbi");
    } else {
       regions_ = new char *[nregs + 1]; // +1 to account for NULL
    int i = 0;
    while (i < nregs)
    {
        regions_[i] = static_cast<char*>(malloc(strlen(regs[i]) + 1));
        strcpy(regions_[i], regs[i]);
        i++;
    }
    regions_[i] = NULL;
    }

    itr_ = tbx_itr_querys(tbx_, regions_[current_reg_]);// calls hts_itr_querys ret an itr or NULL if error
}



// DESTRUCTOIR : deallocating some attributes, probably overkill in C++
htsVCF::~htsVCF()
{
   if (regions_ != nullptr) {
     /* This part is causing problems when regions_ made by me
      int i = 0;
      while ( i < nregs_ ) {
        free(regions_[i]); // freeing mallocated strings
        i++; } */ 
    delete regions_;
   }
    tbx_itr_destroy(itr_);
    tbx_destroy(tbx_);
    if (str_.s) free(str_.s);
}

// GETTERS
char * htsVCF::fname() const { return fp_->fn; }
char * htsVCF::line() { return str_.s; }
int htsVCF::nregs() const { return nregs_; }



//Updating str_.s with the next line of data from the file.
bool htsVCF::next()
{
    if (current_reg_ >= nregs_) return false;
    int ret = 0;
    if (current_reg_ == 0) 
    {
        if ((ret = hts_getline(fp_, '\n', &str_)) >= 0) // 2nd parameter unused, but must be '\n' (or KS_SEP_LINE) (else aborting)
        { // hts_getline updates fp_.lineno
            if ( !(!str_.l || str_.s[0]!= 35) ) {// 35 == conf.meta_char TODO : update comment to be accurate
                return true; // TODO : add why true;
            }
        }
        if (ret < -1) throw std::runtime_error("The reading of the file failed"); //-1 on end-of-file; <= -2 on error
    }

    while ( !itr_ ) { // could be false if region invalid
        if ( ++current_reg_ >= nregs_ )  {
            return false;// no more regions
        }
        else {
            Rcpp::Rcout << "Warning : " << regions_[current_reg_ ] << " seems to be an invalid region. \n";
            itr_ = tbx_itr_querys(tbx_, regions_[current_reg_]);
        }
    }

    ret = tbx_itr_next(fp_, tbx_, itr_, &str_);// == hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx)) @return >= 0 on success, -1 when there is no more data, < -1 on error
    // TODO : see if separate reg utile
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
const char ** htsVCF::list_chroms()
{
    const char **seq = NULL;
    int i, nseq = 0;
    seq = tbx_seqnames(tbx_, &nseq); //very weird func who gives nbr of chroms to nseq and array of chroms to seq
    if (!seq) throw std::runtime_error("Failed to get list of sequence names");
    return seq; // TODO : find a way to free seq ????

}