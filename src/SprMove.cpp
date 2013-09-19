#include <array> 

#include "GibbsProposal.hpp"
#include "AdHocIntegrator.hpp"
#include "SprMove.hpp"
#include "Arithmetics.hpp"
#include "BoundsChecker.hpp"

#define NUM_ITER 3 
// #define ONLY_FIRST  


// TODO make branch stuff more optional (performance)

// TODO constructor instead of extract move info 

void SprMove::applyToTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams) const
{
  applyPath(traln, path, blParams); 
}

void SprMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const
{
  Path anotherPath ; 
  getPathAfterMove(traln, path, anotherPath);
  applyPath(traln, anotherPath, params); 
  path.restoreBranchLengthsPath(traln, params); 
}

 
void SprMove::disorientAtNode(TreeAln &traln, nodeptr p) const
{
  sprDisorientPath(traln,p, path);
}


void SprMove::extractMoveInfo(const TreeAln &traln, std::vector<BranchPlain> description, const std::vector<AbstractParameter*> &params)
{
  sprCreatePath(traln, description.at(0), description.at(1), path, params);
} 


AbstractMove* SprMove::clone() const
{
  return new SprMove; 
}



BranchPlain SprMove::getPruningBranchAfterPrune()const 
{
  return BranchPlain(path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 
}

BranchPlain SprMove::getPruningBranchBeforeOuter() const
{
  return path.at(0); 
}

BranchPlain SprMove::getPruningBranchBeforeInner() const 
{
  return path.at(1); 
}

BranchPlain SprMove::getSubtreeBranchBefore(const TreeAln &traln ) const 
{
  return path.at(0).getThirdBranch(traln, path.at(1)); 
}

BranchPlain SprMove::getSubtreeBranchAfter(const TreeAln &traln ) const 
{
  auto numNodes = path.getNumberOfNodes(); 
  
  auto a =   BranchPlain (path.getNthNodeInPath(1) ,path.getNthNodeInPath(numNodes - 1)); 
  auto b = BranchPlain(path.getNthNodeInPath(1), path.getNthNodeInPath(numNodes-2)); 
  
  return a.getThirdBranch(traln,b); 
}

BranchPlain SprMove::getInsertionBranchBefore() const 
{
  return path.at(path.size() -1 ); 
}

BranchPlain SprMove::getInsertionBranchAfterOuter()  const 
{
  return BranchPlain(path.getNthNodeInPath(1), path.getNthNodeInPath(path.getNumberOfNodes()-1)); 
}


BranchPlain SprMove::getInsertionBranchAfterInner() const 
{
  return BranchPlain(path.getNthNodeInPath(1), path.getNthNodeInPath(path.getNumberOfNodes()-2));
}


BranchPlain SprMove::getEvalBranch(const TreeAln &traln) const
{    
  auto b = path.at(0).getThirdBranch(traln, path.at(1)); 
  auto lastBranch = path.at(path.size()-1); 

  auto bA = BranchPlain(b.getPrimNode(), lastBranch.getPrimNode()), 
    bB = BranchPlain(b.getPrimNode(),lastBranch.getSecNode()); 

  assert(bA.exists(traln) && bB.exists(traln)); 

  auto futureRoot = bA.getThirdBranch(traln, bB ); 
  
  return futureRoot; 
}


#if 0 
/**
   @brief applies the branch length multiplier along the path
   (considering the spr has already been applied to the tree)
 */ 
void SprMove::multiplyBranches(TreeAln &traln, Randomness &rand, double &hastings, PriorBelief &prior, double multiplier, std::vector<AbstractPrior*> brPrs)  const 
{  
  assert(path.size() >= 2); 
  int numBranches = traln.getNumBranches(); 
  assert(numBranches == 1 ); 
  assert(brPrs.size() < 2 ); 
  auto brPr = brPrs[0]; 
 
  int sTNode = path.getNthNodeInPath(1); 
  Branch firstBranch = Branch( path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 

  path.multiplyBranch(traln, rand, firstBranch, multiplier, hastings, prior, brPr); 
  
  /* treat all branches except the first 2 and the last one */
  int ctr = 0; 
  for(nat i = 0; i < path.size() ; ++i)
    {
      if(ctr < 2 )
	continue; 

      path.multiplyBranch(traln, rand, path.at(i), multiplier, hastings, prior, brPr); 
    }

  int lastNode = path.getNthNodeInPath(path.getNumberOfNodes()-1),
    s2LastNode = path.getNthNodeInPath(path.getNumberOfNodes()-2); 
  path.multiplyBranch(traln, rand, Branch(sTNode, lastNode), multiplier, hastings, prior, brPr); 
  path.multiplyBranch(traln, rand, Branch(sTNode, s2LastNode), multiplier, hastings, prior, brPr); 
} 
#endif

/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void SprMove::applyPath(TreeAln &traln, const Path &modifiedPath, 
			const std::vector<AbstractParameter*> &params) const 
{
  assert(modifiedPath.size() > 2 ); 

  /* get the subtree ptr */
  
  auto third =  modifiedPath.at(0).getThirdBranch(traln, modifiedPath.at(1)); 
  nodeptr sTPtr = third.findNodePtr(traln); 
  assert(modifiedPath.at(1).hasNode(sTPtr->number )); 
  
  // finds the two nodeptrs adjacent to the subtree  
  auto prOuterPtr = BranchPlain(modifiedPath.getNthNodeInPath(0) ,modifiedPath.getNthNodeInPath(1)).findNodePtr(traln ),
    prInnerPtr = BranchPlain(modifiedPath.getNthNodeInPath(2), modifiedPath.getNthNodeInPath(1)).findNodePtr(traln ); 

  int lastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1),
    s2LastNode = modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2); 
  nodeptr s2lPtr = BranchPlain(s2LastNode, lastNode).findNodePtr(traln ), // s2l 
    lPtr = BranchPlain(lastNode, s2LastNode).findNodePtr(traln ); // l

  nodeptr toBeInserted = prInnerPtr->back;   
  nodeptr toBeInsertedPath = prOuterPtr->back;  

  assert(toBeInserted->number == sTPtr->number && sTPtr->number == toBeInsertedPath->number); 

  /* prune sTPtr */
  traln.clipNode(prInnerPtr, prOuterPtr); 		 
  traln.setBranch(traln.getBranch(prOuterPtr, params), params); 
  sTPtr->next->back = sTPtr->next->next->back = (nodeptr) NULL; 
  
  // and hook up at the reattachment point, mapping the branches 
  traln.clipNode(toBeInsertedPath, lPtr);
  traln.setBranch(traln.getBranch(lPtr, params), params);

  traln.clipNode(toBeInserted, s2lPtr); 
  traln.setBranch(traln.getBranch(toBeInserted, params), params); 

  assert(BranchPlain(sTPtr->number, s2lPtr->number).exists(traln)); 
  assert(BranchPlain(sTPtr->number, lPtr->number).exists(traln )); 
}



/** 
    @brief Gets the description path after the move has been executed. 

    Can be used for reversal => orientation is inverted 

    Does not contain branch lengths.   
 */ 
void SprMove::getPathAfterMove(const TreeAln &traln, const Path &modifiedPath, Path &resultPath) const 
{  
  int subTreeNode =  modifiedPath.getNthNodeInPath(1);   
  auto lastBranch = modifiedPath.at(modifiedPath.size()-1); 

  resultPath.append(BranchPlain( lastBranch.getSecNode() , subTreeNode ) )  ; 
  resultPath.append(BranchPlain( lastBranch.getPrimNode() , subTreeNode ) )  ; 
  
  assert(modifiedPath.size() > 2); 
  
  // insert the non-descriptive branches in reversed order 
  for(int i = int(modifiedPath.size())-2 ; i > 1 ; --i)
    resultPath.append(modifiedPath.at(i)); 

  // fuse the first two branches 
  auto b1 = modifiedPath.at(0),
    b2  = modifiedPath.at(1); 

  std::set<int> ids; 
  ids.insert(b1.getPrimNode()); 
  ids.insert(b2.getPrimNode()); 
  ids.insert(b1.getSecNode()); 
  ids.insert(b2.getSecNode()); 
  
  assert(ids.find(subTreeNode) != ids.end());   
  ids.erase(subTreeNode); 
  
  assert(ids.size() == 2); 

  auto it = ids.begin(); 
  int a = *it; 
  ++it; 
  int b = *it;   

  resultPath.append(BranchPlain(a,b)); 
  assert(modifiedPath.size() == resultPath.size()); 
}


std::ostream& operator<<(std::ostream &out, const SprMove& rhs)
{
  return out << rhs.path;   
} 




void SprMove::sprCreatePath(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, 
			    Path &pathHere,  const std::vector<AbstractParameter*> &params) const
{
#ifdef EFFICIENT
  // looks very clunky 
  assert(0); 
#endif

  auto chosen = movedInto; // description.at(1); 
  auto prunedTree = mover; // description.at(0); 

  nodeptr p = prunedTree.findNodePtr(traln);
  
  pathHere.clear();   
  pathHere.findPath(traln, p, chosen.findNodePtr(traln));
  pathHere.reverse();   

  auto Tmp = pathHere.at(0); 

  nodeptr aNode = Tmp.findNodePtr(traln); 
  if(traln.isTipNode(aNode))
    aNode = Tmp.getInverted().findNodePtr(traln);
  assert(not traln.isTipNode(aNode)); 
  auto b = pathHere.at(0); 
  while(b.equalsUndirected( pathHere.at(0)) || b.equalsUndirected(BranchPlain(p->number, p->back->number))) 
    {
      aNode = aNode->next; 
      b = BranchPlain(aNode->number, aNode->back->number) ;
    }
  pathHere.reverse();
  pathHere.append(b);  
  pathHere.reverse();

  // TODO inefficient
  pathHere.saveBranchLengthsPath(traln, params); 
}

std::vector<BranchLengths> SprMove::proposeBranches(TreeAln &traln, const std::vector<AbstractParameter*> &params, 
					     LikelihoodEvaluator& eval, double &hastings, Randomness& rand, bool isForward)
{
#if 0 

  bool printVerbose = false; 

  auto result=  std::vector<BranchLength>{};
  // notice: move must have been applied already
  // assert(blParams.size() == 1); 
  // auto param = blParams[0]; 

  tout << MAX_SCI_PRECISION; 

  // pruning point 
  auto branchPrimer =  std::vector<BranchLengths>{};

  if(isForward)
    {
      branchPrimer = 
	{ 
	  traln.getBranch(getSubtreeBranchAfter(traln), params) 
#ifndef ONLY_FIRST
	  , traln.getBranch( getPruningBranchAfterPrune()   ,params)
	  , traln.getBranch( getInsertionBranchAfterInner() , params)
	  , traln.getBranch( getInsertionBranchAfterOuter() , params)
	  , traln.getBranch( getOppositeBranch(traln), params)
#endif
	} ; 
    }
  else 
    {
      branchPrimer = 
	{
	  traln.getBranch( getSubtreeBranchBefore(traln), params)
#ifndef ONLY_FIRST
	  , traln.getBranch( getPruningBranchBeforeOuter(), params)
	  , traln.getBranch( getPruningBranchBeforeInner(), params)
	  , traln.getBranch( getInsertionBranchBefore(), params) 
	  , traln.getBranch( getOppositeBranch(traln) , params )
#endif
	}; 
    } 
  const auto branches = branchPrimer; 

  auto branchNames = std::unordered_map<BranchLengths, std::string>{}; 

  if(isForward)
    {
      branchNames = 
	{
	  make_pair(branches[0], "subtree")
#ifndef ONLY_FIRST
	  ,make_pair(branches[1], "pruning"), 
	  make_pair(branches[2], "insertionInner"), 
	  make_pair(branches[3], "insertionOuter"),
	  make_pair(branches[4], "opposite")
#endif

	}; 
    }
  else 
    {
      branchNames = 
	{
	  make_pair(branches[0], "subtree") 
#ifndef ONLY_FIRST
	  ,make_pair(branches[1], "pruningOuter"),
	  make_pair(branches[2], "pruningInner"),
	  make_pair(branches[3], "insertion"),
	  make_pair(branches[4], "opposite")
#endif
	}; 
    }

  
  // various optima 
  auto optimaMap = std::unordered_map< BranchLengths, double >{}; 
  for(int i = 0; i < NUM_ITER; ++i)
    {
      for( auto branch : branches )
	{
	  auto tmp = traln.getBranch(branch.toPlain(), params); 
	  for(auto &param : params)
	    {
	      double nrd1 = 0; 
	      double nrd2 = 0; 
	      auto result = GibbsProposal::optimiseBranch(traln,branch.toOneLength(param),eval, nrd1, nrd2, 30, param  ); 
	      traln.setBranch(result, param); 
	    }

	  if(optimaMap.find(result) != optimaMap.end())
	    optimaMap.erase(optimaMap.find(result)); 
	  optimaMap[result] = nrd2; 
	}
      if(printVerbose)
	tout << std::endl; 
    }
  if(printVerbose)
    tout << std::endl; 

  
  if(printVerbose)
    tout << "before: " << std::endl; 
  for(auto branch:  branches)
    {
#if 0 
      if(printVerbose)
	tout << branch << "\t" << branch.getInterpretedLength(traln, params) <<  "\t" << branchNames.at(branch) << std::endl; 
#else 
      assert(0); 
#endif
    }
  
  if(printVerbose)
    tout << "after: " << std::endl; 
  for(auto elem : optimaMap)
    {
#if 0 
      if(printVerbose)
	tout << elem.first << "\t" << elem.first.getInterpretedLength(traln, params) << "\t" << branchNames[elem.first] << std::endl; 
#else 
      assert(0); 
#endif
    }
  
  
  if(printVerbose)
    {
      if(isForward)
	tout << "proposing: " << std::endl; 
      else 
	tout << "evaluating: " << std::endl; 
    }
    
  for(auto &branch : branches)
    {
      auto iter = optimaMap.find(branch); 
      auto proposalResult = GibbsProposal::propose(iter->first.getInterpretedLength(traln, param), iter->second, rand); 
      auto tmp = branch; 
      tmp.setConvertedInternalLength(traln, param, proposalResult[0]);
      if( not BoundsChecker::checkBranch(tmp))
	BoundsChecker::correctBranch(tmp); 

      result.push_back(tmp); 

      auto alpha = proposalResult[1]; 
      auto beta = proposalResult[2]; 

      double hastPart = 0; 
      if(isForward)
	{
	  hastPart -= logGammaDensity(tmp.getInterpretedLength(traln,param), alpha, beta) ; 
#if 0 
	  if(printVerbose)
	    tout << tmp << "\t" << tmp.getInterpretedLength(traln, param) << "\t" << branchNames[branch] <<  "\t" << SOME_SCI_PRECISION << hastPart << MAX_SCI_PRECISION   <<std::endl; 
#else 
	  assert(0); 
#endif
	}
      else 
	{
	  hastPart += logGammaDensity(branch.getInterpretedLength(traln, param), alpha, beta); 
#if 0 
	  if(printVerbose)
	    tout << tmp << "\t" << branch.getInterpretedLength(traln, param) << "\t" << branchNames[branch] <<  "\t" << SOME_SCI_PRECISION << hastPart << MAX_SCI_PRECISION   <<std::endl; 
#else 
	  assert(0); 
#endif
	}
      
      hastings += hastPart;       
    }

  // reset 
  for(auto &branch : branches)
    traln.setBranch(branch, params); 

  return result; 
#else 

  // this is a big todo, do not have time right now 
  assert(0); 
  return std::vector<BranchLengths>{};  
#endif
}


void SprMove::sprDisorientPath(TreeAln &traln, nodeptr p, const Path &pathHere) const 
{  
  assert(pathHere.size() > 2) ; 
  
  int first = (int)pathHere.getNthNodeInPath(0),
    last = (int)pathHere.getNthNodeInPath(pathHere.getNumberOfNodes()-1); 

  if(not pathHere.nodeIsOnPath(p->number) || traln.isTipNode(p)
     || p->number == first || p->number == last )
    return; 

  disorientHelper(traln, p);
  
  sprDisorientPath( traln, p->next->back, pathHere); 
  sprDisorientPath( traln, p->next->next->back, pathHere);
}


// #define DO_INTEGRATE 

void SprMove::integrateBranches( TreeAln &traln, const std::vector<AbstractParameter*> blParams, LikelihoodEvaluator &eval, double &hastings ) const 
{
  double prevLnl = traln.getTr()->likelihood;  

  assert(blParams.size() == 1); 
  auto blParam = blParams[0]; 

  double oldHastings = hastings; 

  auto means = std::vector<double>(); 

  auto prunedSubtree = path.at(0).getThirdBranch(traln, path.at(1)); 

  // branches before 
  auto beforeBranches = std::vector<BranchPlain>(); 
  beforeBranches.push_back(path.at(0));
  beforeBranches.push_back(path.at(1));
  beforeBranches.push_back(prunedSubtree); 
  beforeBranches.push_back(path.at(path.size()-1));

  auto resultVec = std::vector<double>(); 
  for(auto b : beforeBranches)
    {
      BranchLength bl = traln.getBranch(b, blParam); 
#ifdef DO_INTEGRATE
      auto samples = ahInt->integrate(bl, traln, 1000, 10);
      auto result = Arithmetics::getMeanAndVar(samples); 
      resultVec.push_back(result.first); 
#else 
      double d1 = 0., d2 = 0. ; 
      eval.evaluate(traln, b, true); 
      auto bCpy = b.toBlDummy(); 
      auto result = GibbsProposal::optimiseBranch(traln, bl, eval, d1,d2, 30, blParam); 
      bCpy.setLength(result); 
      traln.setBranch(bl, blParam); 
      eval.evaluate(traln, b, true); 
      resultVec.push_back(bCpy.getInterpretedLength(traln, blParam)); 
#endif
      // tout << result.first << "\t" ; 

    }

  applyToTree(traln, blParams);

  auto afterBranches = std::vector<BranchPlain>();
  afterBranches.emplace_back(path.getNthNodeInPath(0), path.getNthNodeInPath(2)); 
  afterBranches.emplace_back(path.getNthNodeInPath(1), path.getNthNodeInPath(path.getNumberOfNodes()-2)); 
  afterBranches.emplace_back(path.getNthNodeInPath(1), path.getNthNodeInPath(path.getNumberOfNodes()-1)); 
  afterBranches.emplace_back(prunedSubtree); 
  
  for(auto b : afterBranches)
    {
      BranchLength bl = traln.getBranch(b,blParam); 
#ifdef DO_INTEGRATE
      auto samples = ahInt->integrate(b,traln,1000,10); 
      auto result = Arithmetics::getMeanAndVar(samples); 
      resultVec.push_back(result.first); 
#else 
      double d1 = 0., d2 = 0. ; 
      eval.evaluate(traln, b, true); 
      auto resultLength = GibbsProposal::optimiseBranch(traln, bl , eval, d1,d2, 30, blParam); 
      bl.setLength(resultLength); 
      traln.setBranch(bl, blParam); 
      eval.evaluate(traln, b, true); 
      resultVec.push_back(bl.getInterpretedLength(traln, blParam)); 
#endif
    }

  double newLnl = eval.evaluate(traln, prunedSubtree, true); 

  revertTree(traln, blParams); 
  hastings = oldHastings; 

  auto bip= path.at(0).getThirdBranch(traln, path.at(1)).getBipartition(traln);
  tout << "PARS\t" << getNniDistance() << "\t" << MAX_SCI_PRECISION <<  newLnl - prevLnl << "\t"; 
  for(auto &v : resultVec)
    tout << v << "\t"; 
  tout << bip << std::endl; 
} 




BranchPlain SprMove::getOppositeBranch(const TreeAln &traln ) const 
{
  auto oneNode = path.getNthNodeInPath(path.getNumberOfNodes()-2); 
  auto nodes = traln.getNeighborsOfNode(oneNode); 
  auto otherNode = 0; 
  for(auto node : nodes)
    {
      if(not path.nodeIsOnPath(node) )
	{
	  assert(otherNode == 0); 
	  otherNode = node; 
	}
    }

  return BranchPlain(oneNode, otherNode); 
}
