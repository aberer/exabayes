#ifndef NO_CACHE_POLICY_HPP
#define NO_CACHE_POLICY_HPP

#include "ArrayPolicy.hpp"

class NoCachePolicy : public ArrayPolicy
{
public: 
  NoCachePolicy(const TreeAln &traln ); 

  virtual void imprintPolicy(const TreeAln &traln, ArrayOrientation &arrayOrient) {} 
  virtual void freeMemory()  {} 
  virtual void prepareForEvaluation(TreeAln &traln, BranchPlain virtualRoot, nat models, ArrayOrientation& arrayOrientation ) {} 
  virtual void accountForRejectionPolicy(TreeAln &traln, const std::vector<bool> &partitions, const std::vector<nat>& invalidNodes, ArrayOrientation &arrayOrient); 
  virtual std::unique_ptr<ArrayPolicy> clone() const ; 

}; 



#endif
