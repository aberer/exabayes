#if 0 
#ifndef PRESCORED_DRAWER
#define PRESCORED_DRAWER



#include <unordered_map>

// #include "Branch.hpp"

class BLPairHash; 
class BLPairEqual; 

typedef std::unordered_map<std::pair<Branch,Branch>, double, BLPairHash, BLPairEqual> Branch2parsScores; 
typedef Branch2parsScores Branch2Weight; 

class PrescoredDrawer
{
public: 
  static Branch2parsScores getPrescoredMap(TreeAln &traln ) ;  
  static Branch2Weight getWeights(); 
  
private:  
  
}; 


class BLPairHash
{ 
public: 
  size_t operator()(const std::pair<Branch,Branch> &a ) const 
  {
    // return std::hash<nat>()(b.getPrimNode()) ^ std::hash<nat>()( b.getSecNode()) ; 
    return std::hash<nat>()(a.first.getPrimNode()) 
      ^ std::hash<nat>()(a.second.getPrimNode())
      ^ std::hash<nat>()(a.first.getSecNode())
      ^ std::hash<nat>()(a.second.getSecNode()); 
  }
}; 

class BLPairEqual
{
public:
  bool operator() (const std::pair< Branch, Branch> &b, const std::pair< Branch, Branch> &a) const 
  {
    return ( a.first.equalsUndirected(b.first)  && a.second.equalsUndirected(b.second)) 
      || ( a.first.equalsUndirected(b.second) && a.second.equalsUndirected(b.first) ) ; 
  }
};



#endif

#endif
