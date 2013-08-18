#include "Bipartition.hpp"

#include <iostream>
#include <limits>


// TODO at some point improve this  
nat Bipartition::numBits = sizeof(nat) * 8; 
nat Bipartition::numBitsMinusOne = numBits - 1 ; // not estethic, but efficient 
nat Bipartition::bitShift = 5; 
nat Bipartition::allOne = std::numeric_limits<nat>::max();
std::vector<nat> Bipartition::perBitMask = {
  1u, 
  1u <<  1, 
  1u <<  2,
  1u <<  3,
  1u <<  4,
  1u <<  5,
  1u <<  6,
  1u <<  7,
  1u <<  8,
  1u <<  9,
  1u <<  10,
  1u <<  11,
  1u <<  12,
  1u <<  13,
  1u <<  14,
  1u <<  15,
  1u <<  16,
  1u <<  17,
  1u <<  18,
  1u <<  19,
  1u <<  20,
  1u <<  21,
  1u <<  22,
  1u <<  23,
  1u <<  24,
  1u <<  25,
  1u <<  26,
  1u <<  27,
  1u <<  28,
  1u <<  29,
  1u <<  30,
  1u <<  31
}; 



Bipartition::Bipartition()
  : bip(0)
  , hash(0)
{
  
}


bool Bipartition::operator==(const Bipartition &rhs) const
{
  if( not ( getHash() == rhs.getHash() && bip.size() == rhs.bip.size()) )
    return false; 
  else 
    {
      auto result = true; 
      for(nat i = 0; result && i < bip.size(); ++i)
	result &= ( bip[i] == rhs.bip[i]) ;  
      return result; 
    }      
} 


Bipartition Bipartition::operator| (const Bipartition &rhs) const
{
  auto result = Bipartition{}; 

  nat longestLength = std::max(bip.size() , rhs.bip.size()); 
  nat shorterLength = std::min(bip.size(), rhs.bip.size()); 
  
  result.bip.resize(longestLength, 0u); 
  // std::cout << "resizing to " << longestLength << std::endl; 
  // std::cout << "now: " << result << std::endl; 

  for(nat i = 0; i < shorterLength; ++i)
    result.bip[i] = bip[i] | rhs.bip[i]; 

  if( bip.size() < rhs.bip.size() )
    {
      for(nat i = shorterLength; i < longestLength; ++i)
	result.bip[i] = rhs.bip[i]; 
    }
  else if(rhs.bip.size() < bip.size())
    {
      for(nat i = shorterLength; i < longestLength; ++i)
	result.bip[i] = bip[i]; 
    }
  
  result.hash = hash ^ rhs.hash; 

  return result; 
}



void Bipartition::initializeWithTaxon(nat pos, nat ranNum)
{
  reserve(pos + 1); 
  set(pos); 
  hash = ranNum; 
} 



void Bipartition::reserve(nat num)
{
  ++num ;
  auto bytesNeeded = num / numBits + ( num % numBits > 0 ? 1 : 0  ) ; 
  if( bip.size() * sizeof(nat) < bytesNeeded )
    bip.resize(bytesNeeded , 0); // meh
}

void Bipartition::set(nat pos)
{
  bip[ pos >> bitShift ] |= perBitMask[ pos &  numBitsMinusOne ]; 
}


bool Bipartition::isSet(nat pos) const
{
  return (bip[ pos >> bitShift ] & perBitMask[ pos & numBitsMinusOne ] )   != 0; 
}

void Bipartition::unset(nat pos)
{
  bip[ pos >> bitShift ] &= ( allOne & perBitMask[ pos & numBitsMinusOne ]) ; 
} 


std::ostream& operator<<(std::ostream& out, const Bipartition& rhs)
{
  for(auto &elem : rhs.bip)
    for(nat i = 0; i < Bipartition::numBits; ++i)
      out <<  ( elem & Bipartition::perBitMask[i] ? "1" : "0" ) ; 
  
  return out;   
}


nat Bipartition::count() const 
{
  nat result = 0; 
  for(auto &elem : bip )
    result += __builtin_popcount(elem);

  return result; 
}


void Bipartition::printVerbose(std::ostream &out, const std::vector<std::string> nameMap) const
{
  nat numTax = nameMap.size() -1 ;   
  bool printInverse = count() > numTax / 2;
  
  bool isFirst = true; 
  
  for(nat i = 0; i < numTax; ++i)
    {
      bool isOne =  isSet(i); 
      if( printInverse ^ isOne )
	{
	  out << (isFirst ? "" : "," ) << nameMap[i+1]; 
	  isFirst = false; 
	}
    }
} 



