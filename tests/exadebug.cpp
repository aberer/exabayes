#include <iostream>
#include "Topology.hpp"
#include "NewickParser.hpp"
#include "SPRMove.hpp"

#include <cassert>

#define ASSERT_TRUE assert





int main(int argc, char *argv[])
{
  auto nw = NewickParser<BareTopology>("((1,2),4,(3,(5,6)));");
  nw.initalize();
  auto t2 = nw.getTree();
  return 0;
}
