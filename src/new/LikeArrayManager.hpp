#ifndef LIKEARRAYMANAGER_H
#define LIKEARRAYMANAGER_H

#include <memory>
#include <vector>
#include <unordered_map>
#include <set>

#include "log_double.hpp"


#include "ArrayReservoir.hpp"

#include "Topology.hpp"

#include "Link.hpp"
#include "bitvector.hpp"

#include "ParamContent.hpp"



class CondLikeSubtree
{
  using ParamContentPtr = std::shared_ptr<ParamContent> ;

public:
  CondLikeSubtree(bitvector bip, double *a = nullptr, std::set< ParamContentPtr > pc = {})
    : _bip{bip}
    , _array{a}
    , _paramContents{pc}
  {
  }
  
  CondLikeSubtree(CondLikeSubtree const &rhs ) = default; 
  CondLikeSubtree& operator=( CondLikeSubtree &&rhs)  = default;   
  CondLikeSubtree& operator=( CondLikeSubtree const& rhs)  = default; 
  
  void releaseArray(ArrayReservoir &res )
  {
    res.deallocate(_array);
    _array = nullptr; 
  }
  
private: 
  bitvector _bip;
  double *_array; 
  std::set< ParamContentPtr > _paramContents;
};


/** 
    currently only strives for the basic attributes of likelihoodevaluator. 

    in the future, it should 
    
    * use a defined number of arrays.
    
    * 
    
*/ 

class LikeArrayManager
{
  using ParamContentPtr = std::shared_ptr<ParamContent> ;

public:
  LikeArrayManager( size_t maxArrays , std::shared_ptr<ArrayReservoir> arRes)
    : _arrayReservoir{arRes}
    , _backupArrays{}
    , _curArrays{}
    , _newArrays{}
    , _unusedArrays{}
    , _maxNumArrays{maxArrays}
    , _root{}
  {}

  log_double evaluate(); 
  
  void imprint();
  void reset();
  void makeVanish(bitvector bv); 
  void makeExist(bitvector bv);
  void tellRoot(Topology::iterator it) { _root = *it; }

  // for debug 
  friend std::ostream& operator<<(std::ostream& s,  LikeArrayManager const& c); 

private:                        // METHODS 
  void flushBackup(std::unordered_map<bitvector,CondLikeSubtree> &map) ;

private:                        // ATTRIBUTES 
  std::shared_ptr<ArrayReservoir> _arrayReservoir;
  
  std::unordered_map<bitvector,CondLikeSubtree> _backupArrays; // arrays backed up 
  std::unordered_map<bitvector,CondLikeSubtree> _curArrays;    // existing arrays consistent with the topology
  std::unordered_map<bitvector,CondLikeSubtree> _newArrays; // arrays computed this round 
  std::unordered_map<bitvector,CondLikeSubtree> _unusedArrays; // have been computed at some time, but not used in the current tree 
  
  size_t _maxNumArrays;
  Link _root; 
};


#endif /* LIKEARRAYMANAGER_H */
