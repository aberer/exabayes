#include <streambuf>

using namespace std; 

class teebuf : public streambuf
{
public: 
  teebuf(streambuf *_sb1, streambuf *_sb2)
    : sb1(_sb1)
    , sb2(_sb2)
  {}

private: 
  virtual int overflow(int c)
  {
    if( c == EOF)      
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
    int const r1 = sb1->pubsync();
    int const r2 = sb2->pubsync();
    return r1 == 0 && r2 == 0 ? 0 : -1; 
  }

  streambuf *sb1; 
  streambuf *sb2; 
}; 
