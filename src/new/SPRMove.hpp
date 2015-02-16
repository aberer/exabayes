#ifndef SPRMOVE_H
#define SPRMOVE_H

#include <cassert>
#include <tuple>

#include "extensions.hpp"
#include "Move.hpp"

#include "BareTopology.hpp"




class SPRMove : public Move 
{
  iterator _movedSubtree;
  iterator _regraftLocation;

  iterator _n1;
  iterator _n2; 

public:
  SPRMove( iterator movedSubtree, iterator regraftLocation   ); 
  
  void initialize( )   ; 
  virtual ~SPRMove(){}
  virtual void apply( ) ; 
  
  virtual std::unique_ptr<Move> clone() const; 

  virtual std::unique_ptr<Move> getInverse() const;
  virtual iterator getReferenceBranch()  const; 


  iterator getMovedSubtree() const {  return _movedSubtree; }
  iterator getRegraftLocation() const { return _regraftLocation; }
  
  // must always be the first! 
  Link getBranchAfterPruning() const { return _emergingLinks[0]; }

private:
  void clear() ; 
};



#endif /* SPRMOVE_H */
