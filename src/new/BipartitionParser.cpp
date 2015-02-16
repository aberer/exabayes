#include "BipartitionParser.hpp"

#include <cassert>
#include <algorithm>

#include <iostream>


using std::make_pair; 
using std::pair;
using std::string; 
using std::vector; 



void BipartitionParser::parsePart( vector<anot_tok_t>::const_iterator begin, vector<anot_tok_t>::const_iterator end)
{
  
}

BipartitionParser::v_iter BipartitionParser::findMatchingBracket( v_iter begin, v_iter totalEnd ) const 
{
  if(begin->first != NewickToken::OPEN_BRACKET)
    throw IllegalNewickFormatException("unexpected symbol"); 

  auto ctr = 0;

  while( begin != totalEnd)
    {
      if( begin->first == NewickToken::OPEN_BRACKET)
        ++ctr;
      else if( begin->first == NewickToken::CLOSING_BRACKET)
        --ctr;
      
      if(ctr == 0 )
        break;

      ++begin; 
    }

  if( ctr > 0 )
    throw IllegalNewickFormatException("more opening brackets than closing brackets");

  // go to next element, that is not the bracket   
  ++begin;
  
  auto next = begin + 1; 
  if( next != totalEnd && next->first == NewickToken::EXT_INFO )
    begin = next; 
  
  return begin; 
}




vector<pair<bitvector,string>> BipartitionParser::parseBip(v_iter begin, v_iter end) const 
{
  auto result = vector<pair<bitvector,string>>{};
  
  if(begin->first == NewickToken::NODE_LABEL)
    {
      auto tax = std::stoi(begin->second);
      auto bip = bitvector{tax-1};
      bip.resize(_numTax);

      auto nextOne = begin + 1 ;
      auto str = (nextOne != end && nextOne->first == NewickToken::EXT_INFO)
        ? nextOne->second
        : "";

      result.push_back( make_pair( bip, str) );
    }
  else
    {
      if( begin->first !=   NewickToken::OPEN_BRACKET)
        throw new IllegalNewickFormatException("expected open bracket at "  + begin->second);

      // begin is now inside the brackets, we have to process all elements in here 
      ++begin; 

      auto myBip = bitvector();
      myBip.resize(_numTax);
      
      while( begin != end )
        {
          bool addElems = true; 
          auto partResult = vector<pair<bitvector,string>>{};
          
          if( begin->first == NewickToken::OPEN_BRACKET)
            {
              auto endHere =  findMatchingBracket(begin, end);
              partResult = parseBip(begin, endHere);
              begin = endHere;  // advance
              // std::cout << "end is " << endHere->second << std::endl;
            }
          else if(begin->first == NewickToken::NODE_LABEL)
            {
              auto tmpEnd = begin + 1 ;
              if( tmpEnd != end && tmpEnd->first == NewickToken::EXT_INFO )
                ++tmpEnd;
              partResult = parseBip(begin,tmpEnd); 
              begin = tmpEnd ;  // advance
            }
          else if( begin->first == NewickToken::CLOSING_BRACKET)
            {
              // this means we are done
              auto next = begin + 1;
              auto str = string("") ;
              
              if(next != end && next->first == NewickToken::EXT_INFO)
                {
                  str = next->second;
                  ++begin;      // skip the info 
                }

              // skip the brackets 
              ++begin; 

              result.push_back( make_pair(myBip, str));
              assert(begin == end); 
              addElems = false; 
            }
          else
            {
              throw IllegalNewickFormatException("problem at "  + begin->second);
            }

          
          if(  addElems )
            {
              // last element in vector was the most recent element
              // processed, we add this to the bipartition that
              // respresents this function call
              // std::cout << "at bip " << myBip  << std::endl;
              myBip |=  partResult.back().first; 
          
              // always add the partial result to the accumulated result
              result.insert(result.end(), partResult.begin(), partResult.end() );
            }
        }
    }

  return result; 
}


vector<pair<bitvector,string>>  BipartitionParser::parse(std::vector< anot_tok_t > const& input) 
{
  _numTax = std::count_if(begin(input), end(input),
                          []( pair<NewickToken, string> const& elem){ return elem.first == NewickToken::NODE_LABEL; });
  return  parseBip( input.cbegin(), input.cend() );
}
