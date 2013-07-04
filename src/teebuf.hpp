#ifndef _TEE_BUF_H
#define _TEE_BUF_H

#include <streambuf>


class teebuf : public std::streambuf
{
public: 
  teebuf(std::streambuf *_sb1, std::streambuf *_sb2)
    : sb1(_sb1)
    , sb2(_sb2)
    , isDisabled(false)
  {}

  void disable() {isDisabled = true; }
  

private: 
  virtual int overflow(int c)
  {
    if(isDisabled)
      return EOF; 
    if(  c == EOF)      
      return !EOF; 
    else 
      {	
	int const r1 = sb1->sputc(c); 
	int const r2 = sb2->sputc(c);
	return r1 == EOF || r2 == EOF ? EOF : c ; 
      }
  }

  virtual int sync()
  {
    if(not isDisabled)
      {
	int const r1 = sb1->pubsync();
	int const r2 = sb2->pubsync();
	return r1 == 0 && r2 == 0 ? 0 : -1; 
      }
    return 0; 
  }

  std::streambuf *sb1; 
  std::streambuf *sb2; 
  bool isDisabled; 
}; 


#endif
