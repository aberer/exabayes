/** 
    @brief a helper object making the code concerning a SPR move
    re-usable.
    
    @todo still could use some cleanup.
 */ 

#ifndef _SPR_MOVE_H
#define _SPR_MOVE_H

#include "axml.h"
#include "PriorBelief.hpp"
#include "Path.hpp"


class SprMove
{
public:  
  void destroyOrientationAlongPath( Path& path, tree *tr,  nodeptr p); 
  void multiplyAlongBranchESPR(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior,  Path &modifiedPath, double multiplier, shared_ptr<AbstractPrior> brPr); 
  void applyPathAsESPR(TreeAln &traln, Path &modifiedPath ); 
  void resetAlongPathForESPR(TreeAln &traln, PriorBelief &prior , Path &modifiedPath);   
  void getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath); 
}; 


#endif
