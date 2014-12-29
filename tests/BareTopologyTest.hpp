#include "BareTopology.hpp"
#include "bitvector.hpp"
#include <iostream>


class BareTopologyTest : public ::testing::Test
{
public:
  BareTopologyTest()
    : t{}
  {}

protected:
  virtual void SetUp()
  {
    t = BareTopology();
    auto iter = t.insert(t.begin());
    iter = t.insert(t.begin()); 
    iter = t.insert(t.begin());
    iter = t.begin().next().next().next();

    auto p1 = t.begin();
    auto p2 = t.begin() + 1 ;
    auto p3 = t.begin() + 2 ;

    t.insert(p1);
    t.insert(p2);
    t.insert(p3);
  }

  virtual void TearDown()
  {
  }

protected:
  BareTopology t; 
};


TEST_F(BareTopologyTest, base )
{
  // test if we traverse all taxa, with different starting points 
  for(auto i = 1 ; i < 7 ; ++i )
    {
      bitvector bv ; 
      for(auto iter = t.begin(i); iter != t.end(); ++iter)
        {
          if( iter->isOuterBranch() )
            {

              auto val = iter->getTaxonNode();
              bv.set( val );
            }
        }
      ASSERT_EQ(bv.count(), 6); 
    }

  bitvector bv ; 
  for(auto iter = t.begin(); iter != t.end(); ++iter)
    {
      if( iter->isOuterBranch() )
        {
          
          auto val = iter->getTaxonNode();
          bv.set( val );
        }
    }
  ASSERT_EQ(bv.count(), 6); 
}

