#ifndef NEWICKPARSERTEST_H
#define NEWICKPARSERTEST_H


#include "NewickParser.hpp"


class NewickParserTest : public testing::Test
{
public:
  NewickParserTest()
    : _nw("(1,(2,4),(3,8));")
  {}
  virtual ~NewickParserTest(){}

  virtual void SetUp()
  {
    _nw.tokenizeInput();
    // _nw.printTokens(); 
  }

protected:
  NewickParser<Topology> _nw; 
};
 


TEST_F(NewickParserTest, base)
{
}



#endif /* NEWICKPARSERTEST_H */


