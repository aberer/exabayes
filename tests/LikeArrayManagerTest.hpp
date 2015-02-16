#include "LikeArrayManager.hpp"
#include "NewickParser.hpp"
#include "ArrayReservoir.hpp"
#include "SPRMove.hpp"


TEST(LikeArrayManagerTest, base)
{
  auto nw = NewickParser<Topology>("(1,(2,((8,7),(9,10))),((3,4),(5,6)));");
  nw.initalize();
  auto t = nw.getTree();

  auto res = std::make_shared<ArrayReservoir>(false);
  auto lam = LikeArrayManager( 2 * t.innerSize(), res);

  t.setObservingLnlArrayManager(&lam);

  std::cout << lam << std::endl;
    
  auto move = SPRMove(t.begin(1).opposite(),t.begin(10).opposite().neighbor());
  auto revMove = t.move(move);

  std::cout << lam << std::endl;
}
