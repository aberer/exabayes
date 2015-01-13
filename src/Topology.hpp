#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "BareTopology.hpp"

#include "bitvector.hpp"

#include <set>
#include <unordered_map>


class Topology : public BareTopology 
{
  friend class iterator;

  friend BareTopology::iterator::difference_type distance (iterator first, iterator last); 
  
public:
  Topology();
  
  Topology( std::vector<bitvector> const& bipartitions )
    : BareTopology()
    , _bvs{}
  {
    this->initializeWithBipartitions(bipartitions);
  }

  bool isEquivalent( Topology const& rhs) const; 

  virtual ~Topology() {}
  virtual iterator insert(iterator it, node_id givenId = 0 );
  virtual iterator erase(iterator it);
  virtual std::unique_ptr<Move> move(Move &move ); 

  friend std::ostream& operator<<(std::ostream& s, const Topology& c); 

  bool verifyBipartitions() const ;

private:
  void checkedInsert( Link link, bitvector bip);

  void reorient_helper(iterator it, bool forceExistence);
  void reorient(iterator it, bool forceExistence);

  bitvector getBipOrDummy(iterator it) const; 

private:
  std::unordered_map<Link,bitvector> _bvs; 
};

#endif /* TOPOLOGY_H */
