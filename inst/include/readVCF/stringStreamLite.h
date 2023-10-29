#include <string>

#ifndef _STRINGSTREAMLITE_
#define _STRINGSTREAMLITE_

// to stream strings delimited by single characters.
// the string is modified in-place, make a copy if needed !
// !! don't deallocate the string !!
class stringStreamLite {
  private:

  char * debut;
  char delim;
  char * token;
  bool eof;
  int token_length;

  bool restore;

  public:
  stringStreamLite(char * debut_, char delim_) : debut(debut_), delim(delim_), eof(*debut == 0), restore(false) {}

  stringStreamLite(std::string & s, char delim_) : debut(&s[0]), delim(delim_), eof(*debut == 0), restore(false) {}

  // insère un 0 à la fin du token, renvoie sa longueur,
  // positionne a au début du token suivant (sauf si fin de chaine)
  int next_token() {
    if(restore) { 
      // on remet le délimiteur qui a été remplacé par un zero à l'appel précédent
      debut[-1] = delim;
    }
    token = debut;
    if(*debut == 0) {
      eof = true;
      return 0; 
    }
    // on déplace le début jusqu'à trouver le délimiteur ou 0
    while(*debut != delim && *debut != 0) {
      debut++;
    }

    token_length = (debut-token);

    if(*debut == delim) {
      *debut = 0;
      debut++;
      restore = true;
    } else {
      restore = false;
    }

    return token_length;
  }

  stringStreamLite & operator>>(std::string & s) {
    token_length = next_token();
    s.assign(token);
    return *this;
  }

  stringStreamLite & operator>>(int & x) {
    token_length = next_token();
    if(token_length > 0) {
      x = atoi(token);
    } else {
      x = -2147483648; // R NA_integer_
    }
    return *this;
  } 

  stringStreamLite & operator>>(double & x) {
    token_length = next_token();
    if(token_length > 0) {
      x = atof(token);
    } else {
      x = NAN;
    }
    return *this;
  } 

  operator bool() { 
    return !eof;
  }
};

#endif
