#ifndef MOVE_H
#define MOVE_H

#include <vector>
#include <memory>

#include "Link.hpp"

#include "BareTopology.hpp"


class Move
{
public: 
  using iterator = BareTopology::iterator; 
  
protected: 
  bool _isApplied;

  std::vector<Link> _emergingLinks ;
  std::vector<Link> _vanishingLinks;
  
public:
  Move()
    : _isApplied{false}
    , _emergingLinks{}
    , _vanishingLinks{}
  { }
  
  virtual ~Move(){}

  virtual void apply( ) =0  ; 
  virtual std::unique_ptr<Move> clone() const = 0 ;
  virtual std::unique_ptr<Move> getInverse() const = 0;

  /** 
      @brief gets a many-purpose branch 
      @todo still a bit nebulous, could be the eval branch  
   */ 
  virtual iterator getReferenceBranch()  const = 0; 
  
  std::vector<Link> getEmergingLinks () const { return _emergingLinks; }
  std::vector<Link> getVanishingLinks() const { return _vanishingLinks; }
};



#endif /* MOVE_H */
