#include <Rcpp.h>
#include "htsVCF.h"

// TODO : add other functions called in c°
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
    // TODO : à supp 
    args_t args_hardcode;
    args_hardcode.regions_fname = 0x0;
    args_hardcode.targets_fname = 0x0;
    args_hardcode.print_header = 0;
    args_hardcode.header_only = 0; // TODO : this is what should be changed for f° samples
    args_hardcode.cache_megs = 10; 
    args_hardcode.download_index = 1;
    args_hardcode.separate_regs = 0;
    // TODO : have a look at that
    tbx_conf_t conf_hardcode_instance = { 2, 1, 2, 0, 35, 0 };
    tbx_conf_t *conf_hardcode = &conf_hardcode_instance;

    //originally protected in an if(args_hardcode->cache_megs)
    if (fp_->format.compression == bgzf)
        // bgzf_set_cache_size(hts_get_bgzfp(fp), args_hardcode.cache_megs * 1048576);// TODO : check what hts_get_bgzfp does

    tbx_ = tbx_index_load3(fname, NULL, args_hardcode.download_index ? HTS_IDX_SAVE_REMOTE : 0);//TODO : check to simplify
    if (!regs || !(*regs))// TODO : Peut-être un ternaire pour give regions qd même ? 
    {
        Rcpp::Rcout << "No region specified so will be going through the entire file\n If you wish filter it, please enter a region like so : Chrom_id:First_Index-Last_Index \n";
        //getting an array of all chroms presents in the cvf file, and doing as if these were all the regions asked
        //so putting them in regions_ and updating nregs at the same time
        regions_ = const_cast<char **>(tbx_seqnames(tbx_, &nregs_));//gives back a const char**, but regions_ not const
        if (!regions_)
        fprintf(stderr, "Couldn't get list of sequence names");
    }

    itr_ = tbx_itr_querys(tbx_, regions_[current_reg_]);// calls hts_itr_querys ret an itr or NULL if error

}



// DESTRUCTOIR : deallocating some attributes, probably overkill in C++
htsVCF::~htsVCF()
{
    tbx_itr_destroy(itr_);
    tbx_destroy(tbx_);
    // TODO : put a static variable to see if str was initialised ? 
    //free(str_->s);
}

//un getter de fname,
char * htsVCF::fname() const { return fp_->fn; }//TODO : test if working and accurate

//un getter de line, return str.s
char * htsVCF::line() { return str_.s; }//TODO : test if working 

//un getter de nregs_ to debug
int htsVCF::nregs() const { return nregs_; }



//fonction next qui retourne un booléen et update un membre line
bool htsVCF::next()
{
    
    if (current_reg_ >= nregs_) return false;
    int ret, found = 0; // TODO : wtf

    while ( !itr_ ) { // could be false if region invalid
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
    // if (ret > -1) //throwerror see dropcardile ? means reading failed  
    return true;

}

//fonction qui ne renvoie que les chromosomes présent dans le file
const char ** htsVCF::list_chroms()
{
    const char **seq = NULL;
    int i, nseq = 0;
    seq = tbx_seqnames(tbx_, &nseq); //very weird func who gives nbr of chroms to nseq and array of chroms to seq
    if (!seq) fprintf(stderr, "Couldn't get list of sequence names");
    for (i=0; i<nseq; i++) {
        if (printf(" This is a chrom in your file %s\n", seq[i]) < 0) fprintf(stderr, "Failed to write to stdout");//TODO : rm when tested printing them as a bonus, can be removed later on
    }
    return seq;// TODO : where to put the free(seq) ? Needed if no malloc ? NOPE so all good

}

/* fonction print_header probablement la plus chiante, revoir si je la réécris ou si je reregarde dans htslib ce qu'il fait
modify args check how it works in final */

//check ses fonctions parsing, voir si intégrables ? 
