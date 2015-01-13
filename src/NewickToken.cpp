#include "NewickToken.hpp"
#include <cassert>


std::string tok2str(NewickToken tok )
{
  switch(tok )
    {
    case NewickToken::NODE_LABEL:
      return "NODE_LABEL"; 
    case NewickToken::BRANCH_LENGTH:
      return "BRANCH_LENGTH" ; 
    case NewickToken::EXT_INFO:
      return "EXT_INFO" ; 
    case NewickToken::OPEN_BRACKET:
      return "OPEN_BRACKET"; 
    case NewickToken::CLOSING_BRACKET :
      return "CLOSING_BRACKET" ;
    default:
      assert(0);
      return ""; 
    }
}
