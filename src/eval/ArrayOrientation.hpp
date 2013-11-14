#ifndef _ARRAY_ORIENTATION_HPP
#define _ARRAY_ORIENTATION_HPP

#include <vector>
#include <iostream>

#include "TreeAln.hpp"

#define INVALID 0 

class ArrayOrientation
{
public: 
  ArrayOrientation(const TreeAln &traln); 

  bool isCorrect(nat part, nat id, nat value) const  { return orientation.at(part).at(id) ==  value; }

  nat getOrientation(nat part, nat id ) const ; 

  void setOrientation(nat part, nat id, nat value)  { orientation[part][id] = value; }
  void setOrientationForAllPartitions(nat id ,nat value)  ; 

  void setPartitionInvalid(nat part)  ; 
  void setInvalid(nat part, nat id); 
  
  // void extractOrientation(const TreeAln& traln, nat part); 
  // void extractOrientation( const TreeAln &traln ); 

  friend std::ostream& operator<<(std::ostream& out, const ArrayOrientation &rhs); 

private: 
  std::vector<std::vector<nat> > orientation ; 
}; 

#endif
