#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "BareTopology.hpp"

#include "bitvector.hpp"

#include <set>

#include <unordered_map>

using std::set;
using std::unordered_map; 

using iterator = BareTopology::iterator; 

class Topology : public BareTopology 
{
  friend class iterator;
  
public:
  Topology(); 
  virtual ~Topology() {}
  virtual iterator insert(iterator it);

  friend std::ostream& operator<<(std::ostream& s, const Topology& c); 

  bool verifyBipartitions() const ; 
  
private:
  void checkedInsert( Link link, bitvector bip);

  void reorient_helper(iterator it);
  void reorient(iterator it);

  bitvector getBipOrDummy(iterator it) const; 
    

private:
  unordered_map<Link,bitvector> _bvs; 
};



#endif /* TOPOLOGY_H */
