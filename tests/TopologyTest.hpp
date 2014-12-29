#ifndef TOPOLOGYTEST_H
#define TOPOLOGYTEST_H

#include "Topology.hpp"


TEST(TopologyTest, base )
{
  auto t = Topology();
  t.insert(t.begin());
  t.insert(t.begin());
  t.insert(t.begin());

  t.verifyBipartitions();
  std::cout << t << std::endl;

  
  t.insert(t.begin(1));
  t.verifyBipartitions();
  std::cout << t << std::endl;
    
  t.insert(t.begin(2));
  t.verifyBipartitions();
  std::cout << t << std::endl;

  t.insert(t.begin(3));
  t.verifyBipartitions();
  std::cout << t << std::endl;

  t.insert(t.begin().neighbor() );
  t.verifyBipartitions();
  std::cout << t << std::endl;

  // for(auto pos = t.begin(); pos != t.end(); ++pos)
  //   {
      
  //   }
}





#endif /* TOPOLOGYTEST_H */

