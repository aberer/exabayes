#ifndef NEWICKTOKEN_H
#define NEWICKTOKEN_H

#include <string>

enum class NewickToken
{
  NODE_LABEL,
    BRANCH_LENGTH,
    EXT_INFO,
    OPEN_BRACKET,
    CLOSING_BRACKET
}; 

std::string tok2str(NewickToken tok ) ; 


class IllegalNewickFormatException : public std::exception
{
public: 
  IllegalNewickFormatException(std::string str)
    : _str{str}
  {}

  virtual const char* what() const noexcept
  {
    return _str.c_str();
  }
  
  private:
  std::string _str; 
}; 



#endif /* NEWICKTOKEN_H */
