#ifndef _VCFsnpInfo_
#define _VCFsnpInfo_

// #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
template<typename chrT>
class VCFsnpInfo {
  public:
  chrT chr;
  int pos;
  std::string id;
  std::string ref;
  std::string alt;
  std::string qual;
  std::string filter;
  std::string info;
};

#endif
