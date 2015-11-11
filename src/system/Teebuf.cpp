#include "Teebuf.hpp" 


Teebuf::Teebuf(std::streambuf *_sb1, std::streambuf *_sb2, std::thread::id masterThread)
  : sb1(_sb1)
  , sb2(_sb2)
  , isDisabled(false)
  , _masterThread(masterThread)
{
}


int Teebuf::overflow(int c)
{
  if( MY_TID != _masterThread)
    return 0; 
  
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



int Teebuf::sync()
{
  if(not isDisabled ) 
    {
      int const r1 = sb1->pubsync();
      int const r2 = sb2->pubsync();
      return r1 == 0 && r2 == 0 ? 0 : -1; 
    }
  return 0; 
}
