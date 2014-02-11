#include "config/ExaBlock.hpp"


double ExaBlock::parseScientificDouble(NxsToken& token) 
{
  double result =  0.0 ;
  auto str = token.GetToken(); 
  token.GetNextToken();
  
  bool sawE =  ( str.back()  == 'e' || str.back() == 'E') ; 
  
  if(token.GetToken().EqualsCaseInsensitive("e"))
    {
      sawE = true; 
      str += token.GetToken(); 
      token.GetNextToken();
    }
  
  if(sawE &&  ( token.GetToken().EqualsCaseInsensitive("+") || token.GetToken().EqualsCaseInsensitive("-") ))
    {
      str += token.GetToken(); 
      token.GetNextToken();
    }
  
  if(sawE)
    {
      str += token.GetToken(); 
      token.GetNextToken(); 
    }

  auto &&iss = std::istringstream{str}; 
  iss >> result ; 
  return result; 
} 







