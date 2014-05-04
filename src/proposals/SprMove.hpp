/** 
    @brief a helper object making the code concerning a SPR move
    re-usable.
 */ 


#ifndef _SPR_MOVE_H
#define _SPR_MOVE_H

#include "data-struct/Path.hpp"
#include <unordered_map>

class PriorBelief; 
class LikelihoodEvaluator; 


typedef std::unordered_map<BranchPlain,std::pair<std::string,nat>> branch2PairNameNum; 

class SprMove
{
public:    		
  void applyToTree(TreeAln &traln, const std::vector<AbstractParameter*> &blParams) const; 
  void applyToTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams, LikelihoodEvaluator& eval, bool considerOuter) const; 
  void revertTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams ) const; 
  void revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator& eval, bool considerOuter) const; 

  
  /**
     @brief get the branch we insert into  
   */ 
  BranchPlain getInsertBranch() const
  {
    return path.at(path.size()-1);
  }
  
  // TODO this tuple is a bad idea =(  => make it 2 values
  void extractMoveInfo(const TreeAln &traln, std::tuple<BranchPlain,BranchPlain> description, const std::vector<AbstractParameter*> &params); 
  void extractBranchesOnly(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, Path &pathHere) const ; 
  BranchPlain getEvalBranch(const TreeAln &traln) const; 
  auto moveBranchProposal(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator& eval, Randomness& rand, bool proposeOuter, double thresh, bool sequential)
    -> std::tuple<std::vector<BranchLengths>,log_double,double> ; 

  /**
     @brief gets designations for all branche sin the move. 

     @param isNni -- indicates whether the nni mapping is applied. If
     false, the mrBayes-mapping (single mapping) isi implied.
   */
  branch2PairNameNum getNames(const TreeAln &traln, bool isNni ) const ; 
  /** 
      @brief gets the number of nni moves to be executed in order to
      achieve this move
  */ 
  nat getNniDistance() const { return path.size( )- 2 ;  }

  void integrateBranches( TreeAln &traln,  const std::vector<AbstractParameter*> blParam, LikelihoodEvaluator &eval, double &hastings ) const ; 

  std::vector<nat> getDirtyNodes(const TreeAln& traln, bool considerOuter) const ; 

  Path& getPathHandle(){return path; }
  
  BranchPlain mapBranchNniStepsAfter(const BranchPlain &branch)const  ; 
  BranchPlain mapBranchSingleMapAfter(const BranchPlain &branch ) const; 

  void printBothMappings() const ; 

  SprMove getInverseMove(const TreeAln &traln, const std::vector<AbstractParameter*> &params) const ;

  friend std::ostream& operator<<(std::ostream &out, const SprMove& rhs); 

  static std::vector<SprMove> getAllUniqueMoves(const TreeAln& traln, nat dist); 

  void setPath(Path _path){path=  _path; }

protected:			// METHODS
  BranchPlain getEvalBranchFromPath(const TreeAln &traln, const Path &pathHere ) const ; 
  std::tuple<std::vector<BranchLengths>,log_double, double> proposeBranches(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator &eval, Randomness& rand, bool proposeOuter, double thresh, bool forward); 
  auto proposeBranchesSequentially(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator &eval, 
					  Randomness& rand, bool proposeOuter, double thresh, bool forward )
    -> std::tuple<std::vector<BranchLengths>,log_double, double>; 
  std::vector<BranchPlain> getInvolvedBranchesInOrder(TreeAln& traln, bool outer, const Path &relPath); 
  bool returnCommonBranchesAfter(const BranchPlain &branch,  BranchPlain &result) const; 
  void sprCreatePath(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, Path &path, const std::vector<AbstractParameter*> &params ) const;
  void applyPath(TreeAln &traln, const Path &modifiedPath, const std::vector<AbstractParameter*> &params ) const; 
  Path getPathAfterMove( const Path &modifiedPath ) const; 

private: 			// ATTRIBUTES
  Path path; 
}; 

#endif
