#ifndef BIPARTITIONPARSER_H
#define BIPARTITIONPARSER_H

#include "bitvector.hpp"
#include <vector>
#include <string>
#include "NewickToken.hpp"

class BipartitionParser
{
  // definitions 
  using anot_bv_t = std::pair<bitvector, std::string >;
  using anot_tok_t = std::pair<NewickToken,  std::string>; 

  using v_iter = std::vector<anot_tok_t>::const_iterator; 

  // members 
  std::vector< anot_bv_t >  _anotBv;
  size_t _numTax; 
  
public:
  BipartitionParser()
    : _anotBv{}
    , _numTax{0 }
  {}

  virtual ~BipartitionParser(){}
  
  std::vector<std::pair<bitvector,std::string>>  parse(std::vector< anot_tok_t > const& input) ;
  std::vector< anot_bv_t > getResult() const {return _anotBv; }

private:
  void parsePart( v_iter begin, v_iter end);
  v_iter findMatchingBracket( v_iter begin, v_iter totalEnd ) const ;
  std::vector< std::pair<bitvector,std::string >> parseBip(v_iter begin, v_iter end ) const ; 
};




#endif /* BIPARTITIONPARSER_H */
