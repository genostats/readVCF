#include <Rcpp.h>
#include "htsVCF.h"

// TODO: test with 



extern "C"{
    #include "query_regions.h" //BEFORE htsFile used in line 6 
    htsFile *hts_opening(const char *, const char *);
}


// CONSTRUCTOR : opening the file and calculating nber of regions to go through 
htsVCF::htsVCF(const char *fname,  char ** regs, int nregs)
: regions_(regs), current_reg_(0), nregs_(nregs), str_({0,0,0})
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
    if (!regs || !(*regs))// TODO : Peut-être un ternaire pour give regions qd même ? 
    {
        Rcpp::Rcout << "No region specified so will be going through the entire file\n If you wish filter it, please enter a region like so : Chrom_id:First_Index-Last_Index \n";
        //getting an array of all chroms presents in the cvf file, and doing as if these were all the regions asked
        //so putting them in regions_ and updating nregs at the same time
        regions_ = const_cast<char **>(tbx_seqnames(tbx_, &nregs_));//gives back a const char**, but regions_ not const
        if (!regions_)
        throw std::runtime_error("Failed to look up regions in .tbi");
    }

    itr_ = tbx_itr_querys(tbx_, regions_[current_reg_]);// calls hts_itr_querys ret an itr or NULL if error
}



// DESTRUCTOIR : deallocating some attributes, probably overkill in C++
htsVCF::~htsVCF()
{
    tbx_itr_destroy(itr_);// ok car initialisé dans tous les cas
    tbx_destroy(tbx_);
    if (str_.s) free(str_.s);
}

// GETTERS
char * htsVCF::fname() const { return fp_->fn; }
char * htsVCF::line() { return str_.s; }//TODO : test if working 
int htsVCF::nregs() const { return nregs_; }



//Updating str_.s with the next line of data from the file.
bool htsVCF::next()
{
    if (current_reg_ >= nregs_) return false;
    int ret, found = 0;// TODO : maybe give up found, or add its use
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
        Rcpp::Rcout << "Warning : " << regions_[current_reg_ ] << " seems to be an invalid region. \n";
        if ( ++current_reg_ >= nregs_ )  {
            return false;// no more regions
        }
        else {
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

/* fonction print_header probablement la plus chiante, revoir si je la réécris ou si je reregarde dans htslib ce qu'il fait
modify args check how it works in final 
// What do I want to give back ? a vector o strings with every header ? 
const char**

pb, also is it possible to go back to the beggining of the file ????? 
sinon est-ce envisageable de recréer un object pour ne parcourir que ses headers sans toucher à celui accessible par le user ? 
Semble excessivement compliqué...

*/

//check ses fonctions parsing, voir si intégrables ? 
// TODO : faire une fonction réinitialisant le next au début du file ???
// TODO : check leaks with valgrind without R ? 

