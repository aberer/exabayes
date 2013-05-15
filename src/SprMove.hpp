#include "axml.h"
#include "PriorBelief.hpp"

class SprMove
{
public:
  
  void destroyOrientationAlongPath( Path& path, tree *tr,  nodeptr p); 
  void drawPathForESPR( TreeAln& traln, Randomness &rand, double stopProp ); 
  void multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior ); 
  void applyPathAsESPR(TreeAln &traln ); 
  void resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior); 


}; 
