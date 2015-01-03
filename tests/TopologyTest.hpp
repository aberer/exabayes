#ifndef TOPOLOGYTEST_H
#define TOPOLOGYTEST_H

#include "Topology.hpp"
#include <algorithm>


TEST(TopologyTest, base )
{
  auto t = Topology();
  t.insert(t.begin());
  t.insert(t.begin());
  t.insert(t.begin());

  ASSERT_TRUE(t.verifyBipartitions()) ;
  // std::cout << t << std::endl;

  
  t.insert(t.begin(1));
  ASSERT_TRUE(t.verifyBipartitions());
  // std::cout << t << std::endl;
    
  t.insert(t.begin(2));
  ASSERT_TRUE(t.verifyBipartitions());
  // std::cout << t << std::endl;

  t.insert(t.begin(3));
  ASSERT_TRUE(t.verifyBipartitions());
  // std::cout << t << std::endl;

  t.insert(t.begin().neighbor() );
  ASSERT_TRUE(t.verifyBipartitions());
  // std::cout << t << std::endl;
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
  std::generate( seq.begin(), seq.end(), [&](  ) { auto result = *anIt ; ++anIt;  return result; });

  for(auto &v : seq)
    std::cout << "(" << v << ")\t";
  std::cout << std::endl;
  
  for(auto iter = _topo.begin(); iter != _topo.end();   ++iter)
    {
      std::cout << ctr << "\t" << iter.get() << std::endl;      
      auto tmp = _topo.begin() + ctr ;
      if(tmp != iter )
        {
          std::cout <<  "expected " << iter.get() << " got " << tmp.get() << std::endl;
        }

      ASSERT_TRUE(tmp == iter);

      std::cout << "END test, need to jump  "  << siz - ctr << std::endl;
      tmp = tmp + int(siz - ctr ) ;
      
      std::cout << "end is " << tmp.get()   << std::endl;
      ASSERT_TRUE( tmp == _topo.end() );
      ++ctr;
      std::cout << "\n" << std::endl;
    }
}


#endif /* TOPOLOGYTEST_H */

