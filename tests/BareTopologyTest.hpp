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



TEST_F(BareTopologyTest, erase)
{
  auto toRemove = t.begin( int(t.outerSize() )  ).opposite();
  
  auto it = t.erase(toRemove);
  t.insert( it );

  auto tOrig = t; 
  
  it = t.erase( t.begin(int(t.outerSize())).opposite());
  auto tM1 = t; 
  it = t.erase( t.begin(int(t.outerSize())).opposite());
  auto tM2 = t; 
  it = t.erase( t.begin(int(t.outerSize())).opposite());

  t.insert(t.begin(1));

  ASSERT_TRUE(t == tM2);
  t.insert(t.begin(2));
  ASSERT_TRUE(t == tM1); 
  t.insert(t.begin(3));
  ASSERT_TRUE( t == tOrig );

  it = t.erase( t.begin(int(t.outerSize())).opposite());
  it = t.erase( t.begin(int(t.outerSize())).opposite());
  it = t.erase( t.begin(int(t.outerSize())).opposite());
  it = t.erase( t.begin(int(t.outerSize())).opposite());
  it = t.erase( t.begin(int(t.outerSize())).opposite());

  it = t.erase( t.begin(int(t.outerSize())).opposite());
}



TEST_F(BareTopologyTest, distance)
{
  // iniatial test 
  auto iter = t.begin();
  auto iter2 = t.begin() + 3 ;

  ASSERT_EQ(distance( iter, iter2) , 3 );


  // test various conditions 
  for(auto i1 = t.begin(); i1 != t.end(); ++i1)
    {
      auto ctr = 0; 
      for(auto i2 = i1 ; i2 != t.end(); ++i2)
        {
          ASSERT_EQ(  distance(i1, i2), ctr );
          ++ctr; 
        }
    }
}



TEST_F(BareTopologyTest, altInit)
{
  auto t2 =  BareTopology(  std::vector<size_t>{ 0, 1 , 2  } );
  ASSERT_EQ( t2.size(),  t.size() );   
}



TEST_F(BareTopologyTest, decompose)
{
  auto tCpy = t; 
  auto nums = tCpy.decomposeToInsertionNumbers();
  auto t2 = BareTopology(nums);
  ASSERT_TRUE( t == t2);
}



TEST_F(BareTopologyTest, findPath)
{
  std::cout << t << std::endl;

  auto path = t.begin(3).findPathTo(  t.begin(4));
  for(auto &v : path)
    std::cout << v << ",";
  std::cout << std::endl;

  t.dumpConnections();
}
