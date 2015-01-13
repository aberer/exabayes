#ifndef TOPOLOGYTEST_H
#define TOPOLOGYTEST_H

#include "stacktrace.hpp"

#include "SPRMove.hpp"

#include "Topology.hpp"
#include <algorithm>


TEST(TopologyTest, base )
{
  auto t = Topology();
  t.insert(t.begin());
  t.insert(t.begin());
  t.insert(t.begin());

  ASSERT_TRUE(t.verifyBipartitions()) ;

  t.insert(t.begin(1));
  ASSERT_TRUE(t.verifyBipartitions());

  t.insert(t.begin(2));
  ASSERT_TRUE(t.verifyBipartitions());

  t.insert(t.begin(3));
  ASSERT_TRUE(t.verifyBipartitions());

  t.insert(t.begin().neighbor() );
  ASSERT_TRUE(t.verifyBipartitions());
}



class TopologyTestClass : public ::testing::Test
{
public:
  TopologyTestClass()
    : _topo()
  {}
  
  virtual void SetUp()
  {
    _topo = Topology();
    _topo.insert(_topo.begin());
    _topo.insert(_topo.begin());
    _topo.insert(_topo.begin());

    _topo.insert(_topo.begin(1));
    _topo.insert(_topo.begin(2));
    _topo.insert(_topo.begin(3));
  }
  
protected:
  Topology _topo; 
}; 




TEST_F(TopologyTestClass, size)
{
  BareTopology bt =  _topo; 
  auto iter = bt.begin() + bt.size() ;
  ASSERT_TRUE(iter == bt.end());
}


TEST_F(TopologyTestClass, advance)
{
  auto siz = _topo.size();
  auto ctr = 0;

  auto seq = vector<Link>(siz);
  auto anIt = _topo.begin(); 

  for(auto iter = _topo.begin(); iter != _topo.end();   ++iter)
    {
      auto tmp = _topo.begin() + ctr ;
      ASSERT_TRUE(tmp == iter);
      tmp = tmp + int(siz - ctr ) ;
      ASSERT_TRUE( tmp == _topo.end() );
      ++ctr;
    }

  for(auto i = 1u; i <= _topo.outerSize(); ++i)
    {
      ctr = 0; 
      for(auto iter = _topo.begin(i); iter != _topo.end(); ++iter)
        {
          auto tmp = _topo.begin(i) + ctr ;
          ASSERT_TRUE(tmp == iter);

          ++ctr;
        }
    }
}



TEST_F( TopologyTestClass, erase )
{
  _topo.erase(  _topo.begin( static_cast<node_id>(_topo.outerSize())).opposite() );
  ASSERT_TRUE(_topo.verifyBipartitions());
  
  _topo.erase(  _topo.begin(static_cast<node_id>(_topo.outerSize())).opposite() );
  ASSERT_TRUE(_topo.verifyBipartitions());
    
  _topo.erase(  _topo.begin( (node_id) _topo.outerSize()).opposite() );
  ASSERT_TRUE(_topo.verifyBipartitions());

  _topo.erase(  _topo.begin( (node_id) _topo.outerSize()).opposite() );
  ASSERT_TRUE(_topo.verifyBipartitions());
}


TEST_F(TopologyTestClass, distance)
{
  // this is the same test as in BareTopologyTest. However, we test
  // the different bip-based implementation (that uses short-cuts)
  
  // iniatial test 
  auto iter = _topo.begin();
  auto iter2 = _topo.begin() + 3 ;

  ASSERT_EQ(distance( iter, iter2) , 3 );

  // test various conditions 
  for(auto i1 = _topo.begin(); i1 != _topo.end(); ++i1)
    {
      auto ctr = 0; 
      for(auto i2 = i1 ; i2 != _topo.end(); ++i2)
        {
          auto dist =  distance(i1, i2);
          ASSERT_EQ(  ctr  , dist );
          ++ctr;
        }
    }
  
  // now for different starting iterators  
  for(auto i = 1u; i <= _topo.outerSize() ; ++i)
    {
      for(auto i1 = _topo.begin(i); i1 != _topo.end(); ++i1)
        {

          auto ctr = 0; 
          for(auto i2 = i1 ; i2 != _topo.end(); ++i2)
            {
              auto dist =  distance(i1, i2);
              ASSERT_EQ(  dist  , ctr );
              ++ctr; 
            }
        }
    }
}




TEST(TopologyTest, altInsert)
{
  auto t =   Topology();
  t.insert(t.begin(), 4);
  t.insert(t.begin(), 5);
  t.insert(t.begin(), 6);

  t.insert(t.begin(4),1);
  t.insert(t.begin(5),2);
  t.insert(t.begin(6),3);
  
  // std::cout << t << std::endl;
  
  // TODO test for equality 

  // ASSERT_TRUE(false);

  auto t2 = Topology();
  t2.insert(t2.begin());
  t2.insert(t2.begin());
  t2.insert(t2.begin());
  
  t2.insert(t2.begin(1),4);
  t2.insert(t2.begin(2),5);
  t2.insert(t2.begin(3),6);

  // std::cout << t2 << std::endl;

  // ASSERT_TRUE(t == t2 );
  ASSERT_TRUE( t.isEquivalent(t2));
}


#include "NewickParser.hpp"

TEST(TopologyTest, initWithBip)
{
  auto str = "((1,2),(3,4),(5,6));" ;
  // std::cout << "tree: " << str << std::endl;
  auto pars = NewickParser<Topology>( str);
  pars.initalize();
  auto annotBip = pars.getAnnotations();

  auto bips = vector<bitvector>(annotBip.size());
  std::transform( annotBip.begin(), annotBip.end(), bips.begin() ,
                  [](  std::pair<bitvector, std::string > const &elem  ) { return elem.first ; }  );

  // std::cout << "bips are: " << std::endl;
  // for(auto &bip : bips)
    // std::cout << bip << std::endl;

  auto t = Topology (  bips);

  auto t2 = Topology();
  t2.insert(t2.begin());
  t2.insert(t2.begin(),3);
  t2.insert(t2.begin(),5);
  t2.insert(t2.begin(1),2);
  t2.insert(t2.begin(3),4);
  t2.insert(t2.begin(5),6);

  // std::cout << "con t" << std::endl;
  // t.dumpConnections();
  // std::cout << "con t2" << std::endl;
  // t2.dumpConnections();
  
  // std::cout << t2  << std::endl;
  // std::cout << t << std::endl;
  
  ASSERT_TRUE(t.isEquivalent(t2) );
}



TEST_F(TopologyTestClass, move)
{
  auto move1 = SPRMove(_topo.begin(2).opposite(), _topo.begin(1));
  _topo.move( move1 );
  ASSERT_TRUE(_topo.verifyBipartitions());

  auto move2 = ( SPRMove(_topo.begin(5).opposite(), _topo.begin(6) )); 
  _topo.move( move2 );
  ASSERT_TRUE(_topo.verifyBipartitions());

  auto nw = NewickParser<Topology>("((1,2),4,(3,(5,6)));");
  nw.initalize();
  auto t2 = nw.getTree();
  
  ASSERT_TRUE( _topo.isEquivalent(t2) );
}


#endif /* TOPOLOGYTEST_H */
