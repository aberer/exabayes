#ifndef SAMPLEDBIPARTITIONS_H
#define SAMPLEDBIPARTITIONS_H

#include <unordered_map>

#include "Serializable.hpp"
#include "bitvector.hpp"

class ParallelSetup; 
class Topology; 
class TreeAln; 

class SampledBipartitions : public Serializable
{
public:                         // INHERITED 
  virtual void deserialize( std::istream &in )  ; 
  virtual void serialize( std::ostream &out) const;
  
public:
  SampledBipartitions()
    : _bip2trees{}
  {}

  void addTree(TreeAln const& traln, uint64_t trIndex) ;
  void addTree(Topology const &t, uint64_t trIndex); 
  std::unordered_map<bitvector,bitvector> getMap() const { return _bip2trees; }

  void reduce( SampledBipartitions const& others);

  friend std::ostream& operator<<(std::ostream&, const SampledBipartitions& rhs);
  
  
private: 
  std::unordered_map<bitvector,bitvector> _bip2trees;
};


#endif /* SAMPLEDBIPARTITIONS_H */
