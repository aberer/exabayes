#ifndef FULL_CACHE_POLICY_HPP
#define FULL_CACHE_POLICY_HPP

#include "ArrayPolicy.hpp" 
#include "ArrayOrientation.hpp"
#include "ArrayRestorer.hpp"

class FullCachePolicy : public ArrayPolicy
{
public: 
  FullCachePolicy(const TreeAln& traln, bool cacheTipTip, bool cacheTipInner  ); 

  virtual void imprintPolicy(const TreeAln &traln, ArrayOrientation &arrayOrient) ; 
  virtual void freeMemory() ;  
  virtual void accountForRejectionPolicy(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient); 
  virtual void prepareForEvaluation(TreeAln &traln, BranchPlain virtualRoot, nat model  , ArrayOrientation& arrayOrientation ); 

  virtual std::unique_ptr<ArrayPolicy> clone() const ; 
  
  void enableRestoreGapVector() {restorer.enableRestoreGapVector(); } 

private: 
  ArrayRestorer restorer; 
  

}; 

#endif
