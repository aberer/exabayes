#ifndef SPRMOVE_H
#define SPRMOVE_H

#include <tuple>

class SPRMove
{
  iterator _movedSubtree;
  iterator _regraftLocation;

  vector<Link> _emergingLinks ;
  vector<Link> _vanishingLinks;

public:
  SPRMove( iterator movedSubtree, iterator regraftLocation   )
    : _movedSubtree{movedSubtree}
    , _regraftLocation{regraftLocation}
    , _emergingLinks{}
    , _vanishingLinks{}
  {
  }

  
  vector<Link> getEmergingLinks ()
  {
    if( _vanishingLinks.size() == 0 )
      computeVanishEmergeLinks( );
    return _emergingLinks; 
  }

  vector<Link> getVanishingLinks()
  {
    if( _emergingLinks.size() == 0  )
      computeVanishEmergeLinks();
    return _vanishingLinks; 
  }


  Link getBranchAfterPruning() 
  {
    if( _emergingLinks.size() == 0 )
      computeVanishEmergeLinks();

    // must always be the first! 
    return _emergingLinks[0]; 
  }
  
  
private: 
  void computeVanishEmergeLinks( )   
  {
    auto p = _regraftLocation->primary(); 
    auto s = _regraftLocation->secondary();
    
    auto n1 = _movedSubtree.neighbor();
    auto n2 = n1.neighbor();

    auto branchAfterPruning =  Link(n1->secondary() , n2->secondary()); 
    
    _emergingLinks = vector<Link>{ branchAfterPruning , // DONT move this guy 
                                   Link(p, _movedSubtree->primary()),
                                   Link(s, _movedSubtree->primary())};
    
    _vanishingLinks = vector<Link> {*n1, *n2, *_regraftLocation};
  }
  

};



#endif /* SPRMOVE_H */
