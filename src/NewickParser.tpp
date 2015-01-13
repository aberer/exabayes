#include "NewickParser.hpp"
#include "BipartitionParser.hpp"

#include <tuple>
#include <algorithm>
#include <cassert>
#include <exception>

#include <sstream>

using std::vector;
using std::string;

using std::make_pair; 

template<class TREE>
vector<typename NewickParser<TREE>::tok_t> NewickParser<TREE>::tokenizeInput() const
{
  auto result = vector<tok_t>{}; 
  auto iter = _input.begin();
  bool foundEnd = false;
  
  while(iter != _input.end())
    {
      if( isWhiteSpace(*iter ) || *iter == ',')
        {
          ++iter; 
        }
      if( *iter == '(' ) // brackets 
        {
          result.push_back( make_pair(NewickToken::OPEN_BRACKET, string{*iter}) );
          ++iter; 
        }
      else if(  *iter == ')')
        {
          result.push_back( make_pair(NewickToken::CLOSING_BRACKET, string{*iter}) );
          ++iter; 
        }
      else if( isNumeric(*iter )) // a node label  
        {
          auto &&ss = std::stringstream{}; 
          while( iter != _input.end() &&  isNumeric(*iter) )
            {
              ss << *iter ;
              ++iter ; 
            }
          result.push_back(  make_pair(NewickToken::NODE_LABEL, ss.str()));
        }
      else if( *iter == ':')
        {
          ++iter;
          auto &&ss = std::stringstream{};
          while( isPartOfDouble(*iter) )
            ss << *iter ;

          if( not isDouble(ss.str()) )
            throw IllegalNewickFormatException( ">" + ss.str()  + "< is not a valid double" ); 
          
          result.push_back( make_pair(NewickToken::BRANCH_LENGTH, ss.str())) ;
        }
      else if( *iter == '[' )
        {
          ++iter;
          auto &&ss = std::stringstream();
          while( *iter != ']')
            ss << *iter;
          ++iter; 
        }
      else if( *iter == ';')
        {
          ++iter; 
          foundEnd = true;
          while(  iter !=  _input.end() )
            {
              if( not isWhiteSpace(*iter))
                throw IllegalNewickFormatException( "found non-whitespace characters after semicolon delimiter: " + *iter);
              ++iter; 
            }
        }
      else
        {
          std::cout << "unknown character >" <<  *iter << "<" << std::endl;
        }
    }

  if( not foundEnd)
    throw IllegalNewickFormatException("did not find a terminating semicolon");

  return result; 
}


template<class TREE>
bool NewickParser<TREE>::isWhiteSpace( char c) 
{
  // any more? 
  return c == ' ' || c == '\t' || c == '\n'; 
}


template<class TREE>
bool NewickParser<TREE>::isNumeric(char c) 
{
  return  ( char(48) <= c )  && ( c <=  char(57) )  ; 
}

template<class TREE>
bool NewickParser<TREE>::isAlpha(char c) 
{
  return
    ( char(65) <= c && c <= char(90) )
    || ( char(97) <= c &&  c <= char(122) ) ;   
}

template<class TREE>
bool NewickParser<TREE>::isAlphaNumeric(char c)  
{
  return isAlpha(c)  || isNumeric(c); 
}


template<class TREE>
bool NewickParser<TREE>::isScientificDouble( std::string str) 
{
  auto iter = str.begin();
  if(*iter == '-')
    ++iter;
  
  while(isNumeric(*iter))
    ++iter;

  if( *iter == '.')
    ++iter;

  while(isNumeric(*iter))
    ++iter; 

  if( *iter == 'e' || *iter == 'E')
    ++iter;
  if( *iter == '-')
    ++iter;

  while( isNumeric(*iter) )
    ++iter;

  return  iter == str.end();
}


template<class TREE>
bool NewickParser<TREE>::isNormalDouble( std::string str ) 
{
  auto iter = str.begin();
  if(*iter == '-')
    ++iter;
  while( isNumeric(*iter) )
    ++iter;
  if( *iter == '.')
    ++iter;
  while(isNumeric(*iter))
    ++iter;

  return iter == str.end(); 
}


template<class TREE>
bool NewickParser<TREE>::isDouble( std::string str ) 
{
  return isNormalDouble(str) || isScientificDouble(str); 
}


template<class TREE>
bool NewickParser<TREE>::isPartOfDouble(  char c) 
{
  return
    c == 'e'
    || c == 'E'
    || c == '-'
    || c == '.'
    || isNumeric(c); 
}


template<class TREE>
void NewickParser<TREE>::initalize()
{
  auto tokens = tokenizeInput();
  auto bipParser = BipartitionParser ();
  _annotations = bipParser.parse(tokens);
  auto bips = vector<bitvector>();
  std::transform(_annotations.begin(), _annotations.end(), std::back_inserter(bips), [](std::pair<bitvector, std::string > const& elem){ return elem.first; });
  _tree = TREE(bips); 
}


