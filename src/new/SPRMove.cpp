#include "SPRMove.hpp"


SPRMove::SPRMove( iterator movedSubtree, iterator regraftLocation   )
  : Move()
  , _movedSubtree{movedSubtree}
  , _regraftLocation{regraftLocation}
  , _n1{_movedSubtree.neighbor()}
  , _n2{_n1.neighbor() }
{
  assert( movedSubtree.isValid());
  assert( regraftLocation.isValid() );

  // cannot check that here   
  if(  _movedSubtree->isOuterBranch()  && _movedSubtree->primary() == _movedSubtree->getTaxonNode()  )
    throw IllegalTreeOperation("Tried to move entire tree."); 
  
  initialize();
}


void SPRMove::initialize( )   
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


void SPRMove::apply( ) 
{
  _isApplied = true; 
}

  
std::unique_ptr<Move> SPRMove::clone() const
{
  return  make_unique<SPRMove>(_movedSubtree, _regraftLocation); 
}



std::unique_ptr<Move> SPRMove::getInverse() const
{
  assert( _isApplied ) ;
  auto oldInsert  =  iterator(  _movedSubtree.getTopologyHandle(), Link( _n1->secondary() , _n2->secondary() ));
  return make_unique<SPRMove>( _movedSubtree,oldInsert ); 
}


void SPRMove::clear() 
{
  _emergingLinks = std::vector<Link>();
  _vanishingLinks = std::vector<Link>();                           
}


BareTopology::iterator SPRMove::getReferenceBranch()  const
{
  return _movedSubtree; 
}
