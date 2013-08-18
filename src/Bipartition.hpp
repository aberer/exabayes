#ifndef BIPARTITION_H
#define BIPARTITION_H

#include <vector>  
#include <cassert>
#include <ostream>
#include <string>
#include <unordered_map>

#include "common.h"

class Bipartition
{
public: 
  explicit Bipartition();

  void setHash(nat h) {assert(h != 0 ) ; hash = h; }
  nat getHash() const {assert(hash != 0 ) ; return hash; }
  /**
     @brief checks equality of two bipartitions
   */
  bool operator==(const Bipartition &rhs) const; 
  /** 
      @brief performs an or-operation on two bipartitions and returns
      the result
      
      @notice the hash may loose its validity, if the bipartitions are
      or-ed that have an intersection.
   */ 
  Bipartition operator| (const Bipartition &rhs) const;
  /** 
      @brief sets the bit at position pos 
   */ 
  void set(nat pos) ; 
  /** 
      @brief unsets the bit at position pos 
   */ 
  void unset(nat pos); 
  bool isSet(nat pos) const;
  void initializeWithTaxon(nat pos, nat ranNum); 
  /**
     @brief makes the bit vector sufficiently long, s.t. num can be
     set

     @notice do not overuse, this is somewhat expensive
   */ 
  void reserve(nat num);
  /**
     @brief gets the number of bits set
   */ 
  nat count() const ;
  /** 
      @brief prints the bipartition in a readable manner
   */ 
  void printVerbose(std::ostream &outt, const std::vector<std::string> nameMap) const; 

  static std::vector<nat> perBitMask;
  static nat numBits;
  static nat numBitsMinusOne; 
  static nat bitShift; 
  static nat allOne; 

  friend std::ostream& operator<<(std::ostream& out, const Bipartition& rhs); 

private: 
  std::vector<nat> bip; 
  nat hash; 

}; 


namespace std
{
  template<> struct hash<Bipartition>
  {
    size_t operator()(const Bipartition& b) const 
    {
      return b.getHash(); 
    }
  }; 

}

#endif
