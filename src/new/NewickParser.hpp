#ifndef NEWICKTOKENIZER_H
#define NEWICKTOKENIZER_H

#include <string>
#include <vector>
#include <iostream>

#include "bitvector.hpp"
#include "NewickToken.hpp"

template<class TREE>
class NewickParser
{
  friend class NewickParserTest;
  friend class BipartitionParserTest; 
  
  using tok_t = std::pair<NewickToken,  std::string>; 
  using tokenIterator = std::vector< tok_t >::const_iterator;
  using anot_bv_t = std::pair<bitvector, std::string >;
  
public:
  NewickParser()
    : _input()
    , _tree{}
    , _annotations{}
  {}
  
  NewickParser(std::string str)
    : NewickParser()
  {
    _input = str; 
  }

  virtual ~NewickParser(){}
  void initalize();

  TREE getTree() const {return _tree; }
  std::vector<anot_bv_t> getAnnotations() const { return _annotations; }

private:                             // METHODS
  std::vector< tok_t > tokenizeInput() const ;

private:                             // STATIC METHODS
  static bool isNumeric(char c)  ;
  static bool isAlpha(char c)  ;
  static bool isAlphaNumeric(char c)  ;
  static bool isWhiteSpace( char c)   ; 
  static bool isScientificDouble( std::string str)  ;
  static bool isDouble( std::string str ) ;
  static bool isNormalDouble( std::string str ) ; 
  static bool isPartOfDouble(  char c)  ;

private:
  std::string _input;
  TREE _tree;
  std::vector<anot_bv_t> _annotations;
  // this could contain a mapping of taxa to ids in the tree ... 
}; 


#include "NewickParser.tpp"

#endif /* NEWICKTOKENIZER_H */
