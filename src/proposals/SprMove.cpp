#include <array> 

#include <cmath>

#include "GibbsProposal.hpp"
#include "AdHocIntegrator.hpp"
#include "SprMove.hpp"
#include "math/Arithmetics.hpp"
#include "system/BoundsChecker.hpp"

#define NUM_ITER 3 
// #define ONLY_FIRST  

// #define VERBOSE_INFO

// TODO constructor instead of extract move info 

void SprMove::applyToTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams) const
{
  applyPath(traln, path, blParams); 
}

void SprMove::applyToTree(TreeAln &traln,const std::vector<AbstractParameter*> &blParams, LikelihoodEvaluator& eval, bool considerOuter) const
{
  applyToTree(traln, blParams ); 
  auto dirtyNodes = getDirtyNodes(traln, considerOuter); 
  for(auto &elem : dirtyNodes)
    eval.markDirty( traln, elem); 
}



void SprMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params) const
{
  auto anotherPath = getPathAfterMove( path);
  applyPath(traln, anotherPath, params); 
  path.restoreBranchLengthsPath(traln, params); 
}


void SprMove::revertTree(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator& eval, bool considerOuter) const
{
  auto anotherPath = getPathAfterMove( path);
  applyPath(traln, anotherPath, params); 
  path.restoreBranchLengthsPath(traln, params); 
  
  auto dirtyNodes = getDirtyNodes(traln, considerOuter);
  for(auto &elem : dirtyNodes)
    eval.markDirty(traln,elem);      
}


void SprMove::extractMoveInfo(const TreeAln &traln, std::tuple<BranchPlain,BranchPlain> description, const std::vector<AbstractParameter*> &params)
{
  sprCreatePath(traln, std::get<0>(description), std::get<1>(description), path, params);
} 


BranchPlain SprMove::getEvalBranchFromPath(const TreeAln &traln, const Path &pathHere ) const 
{
  auto b = pathHere.at(0).getThirdBranch(traln, pathHere.at(1)); 
  auto lastBranch = pathHere.at(pathHere.size()-1); 

  auto bA = BranchPlain(b.getPrimNode(), lastBranch.getPrimNode()), 
    bB = BranchPlain(b.getPrimNode(),lastBranch.getSecNode()); 

  assert(bA.exists(traln) && bB.exists(traln)); 

  auto futureRoot = bA.getThirdBranch(traln, bB ); 
  
  return futureRoot; 
}



BranchPlain SprMove::getEvalBranch(const TreeAln &traln) const
{    
  return getEvalBranchFromPath(traln, path); 
}


/**
   @brief applies the path onto the tree 

   first two branches in path define subtree, all further branches are
   traversed, last branch is the insertion branch. 
 */
void SprMove::applyPath(TreeAln &traln, const Path &modifiedPath, 
			const std::vector<AbstractParameter*> &params) const 
{
  assert(modifiedPath.size() > 2 ); 

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


Path SprMove::getPathAfterMove( const Path &modifiedPath ) const 
{ 
// create first elem 
  auto resultPath = Path{}; 

  resultPath.append(BranchPlain(modifiedPath.getNthNodeInPath(0) , modifiedPath.getNthNodeInPath(2) ));
  for(nat i = 2; i < modifiedPath.size()-1 ; ++i)
    resultPath.append(modifiedPath.at(i));

  resultPath.append(BranchPlain(modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-2), 
				modifiedPath.getNthNodeInPath(1))); 
  resultPath.append(BranchPlain(modifiedPath.getNthNodeInPath(modifiedPath.getNumberOfNodes()-1), 
				modifiedPath.getNthNodeInPath(1))); 
  resultPath.reverse();

  return resultPath; 
}


std::ostream& operator<<(std::ostream &out, const SprMove& rhs)
{
  return out << rhs.path;   
} 


void SprMove::extractBranchesOnly(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, Path &pathHere) const 
{
  auto chosen = movedInto; 
  auto prunedTree = mover; 

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
}


void SprMove::sprCreatePath(const TreeAln &traln, BranchPlain mover, BranchPlain movedInto, Path &pathHere,  const std::vector<AbstractParameter*> &params) const
{
  extractBranchesOnly(traln, mover, movedInto, pathHere); 
  pathHere.saveBranchLengthsPath(traln, params); 
}


std::vector<BranchPlain> SprMove::getInvolvedBranchesInOrder(TreeAln& traln, bool outer, const Path &relPath)
{
  auto result = std::vector<BranchPlain>{}; 

  if(outer)
    {
      auto b = relPath.at(0).toPlain(); 
      result.push_back(b);
    }

  for(nat i = 1; i < relPath.size() -1; ++i)
    {
      if(outer)
	{
	  auto b = relPath.at(i-1).getThirdBranch(traln, relPath.at(i)).toPlain();
	  result.push_back(b);
	}
      result.push_back(relPath.at(i));
    }
  
  if(outer)
    {
      auto b = relPath.at(relPath.size() - 1 ); 
      result.push_back(b);

      auto b1 = relPath.at(relPath.size()-1); 
      auto b2 = relPath.at(relPath.size()-2).toPlain(); 
      auto c = b1.getThirdBranch( traln, b2);
      result.push_back(c);
    }
  return result; 
}


auto SprMove::moveBranchProposal(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator& eval,
				 Randomness& rand, bool proposeOuter, double thresh, bool sequential)
  -> std::tuple<std::vector<BranchLengths>,log_double, double> 
{
  auto partResult = 
    sequential
    ? proposeBranchesSequentially(traln, params, eval, rand, proposeOuter, thresh, true)
    : proposeBranches(traln, params, eval, rand, proposeOuter, thresh, true);

  for(auto b : std::get<0>(partResult) )
    {
      auto bl = b.toBlDummy(); 
      bl.setLength(b.getLengths()[0]); 
    }

  applyToTree(traln,params, eval, proposeOuter);
  auto inverseMove = getInverseMove(traln, params); 

  auto backPartResult = sequential
    ? inverseMove.proposeBranchesSequentially(traln, params, eval,rand,proposeOuter, thresh, false)
    : inverseMove.proposeBranches(traln, params, eval,rand,proposeOuter, thresh, false); 
  revertTree(traln, params, eval, proposeOuter);

  double maxImpact = std::get<2>(partResult);
  return make_tuple(std::get<0>(partResult), std::get<1>(backPartResult)  / std::get<1>(partResult), maxImpact);
}



auto SprMove::proposeBranchesSequentially(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator &eval, 
					  Randomness& rand, bool proposeOuter, double thresh, bool forward )
  -> std::tuple<std::vector<BranchLengths>,log_double, double>
{
  log_double probability = log_double::fromAbs(1.);
  auto result = std::vector<BranchLengths>{}; 

  assert(params.size() == 1 ); 
  auto param = params[0];
  
  applyToTree(traln, params, eval, proposeOuter); 
  auto backPath = getPathAfterMove(path);
  
  auto involvedBranches = getInvolvedBranchesInOrder(traln, proposeOuter, backPath); 
  
  auto branch2lengthOrig = std::unordered_map<BranchPlain,double>{}; 
  for(auto b : involvedBranches)
    branch2lengthOrig[b.toPlain()] = traln.getBranch(b.toPlain(), param).getLength(); 

  std::reverse(involvedBranches.begin(), involvedBranches.end()); 

  double ratio = 0; 
  for(auto b : involvedBranches)
    {
      auto bl = traln.getBranch(b, param); 

      double before = bl.getInterpretedLength(traln, param);

      auto optTuple =  GibbsProposal::optimiseBranch(traln, bl, eval, 30, param); 
      bl.setLength(optTuple[0]); 
      auto nrd1 = optTuple[1]; 
      auto nrd2 = optTuple[2]; 

      auto proposal = GibbsProposal::propose(bl.getInterpretedLength(traln, param), nrd1, nrd2, rand); 
      auto proposedLength = proposal[0]; 

      double after = proposedLength; 
    
      double impactHere = fabs(log(after) - log(before)); 
      if(ratio < impactHere)
	ratio = impactHere;

      bl.setConvertedInternalLength(traln, param, proposedLength);
      
      if(not BoundsChecker::checkBranch(bl))
	BoundsChecker::correctBranch(bl); 

      auto probPart = Density::lnGamma(bl.getInterpretedLength(traln,param), proposal[1], proposal[2]);       

      traln.setBranch(bl,param); 

      auto bls = bl.toBlsDummy(); 
      bls.setLengths({bl.getLength()});
      result.push_back(bls); 

      probability *= probPart; 
    }

  for(auto elem : branch2lengthOrig)
    {
      auto b = elem.first.toBlDummy(); 
      b.setLength(elem.second); 
      traln.setBranch(b, param);
    }
  revertTree(traln,params, eval, proposeOuter );   
  
  return std::make_tuple(result, probability, ratio);
}


std::vector<BranchPlain> findAdjacentBranches(const BranchPlain& branch, const std::vector<BranchPlain> branches)
{
  auto result = std::vector<BranchPlain>{}; 
  for(auto &b : branches)
    {
      if(not branch.equalsUndirected(b) && b.isAdjacent(branch))
	{
	  result.push_back(b);
	}
    }

  return result; 
}


auto SprMove::proposeBranches(TreeAln &traln, const std::vector<AbstractParameter*> &params, LikelihoodEvaluator &eval, 
			      Randomness& rand, bool proposeOuter, double thresh, bool forward )
  -> std::tuple<std::vector<BranchLengths>,log_double, double>
{
  auto probability = log_double::fromAbs(1.); 
  auto result = std::vector<BranchLengths>{};

  assert(params.size() == 1); 
  auto param = params[0]; 

  applyToTree(traln, params, eval, proposeOuter); 

  auto backPath = getPathAfterMove(path);

  auto involvedBranches = getInvolvedBranchesInOrder(traln, proposeOuter, backPath); 

  std::reverse(involvedBranches.begin(), involvedBranches.end()); 
  
  auto branch2lengthOrig = std::unordered_map<BranchPlain,double>{}; 
  for(auto b : involvedBranches)
    branch2lengthOrig[b.toPlain()] = traln.getBranch(b.toPlain(), param).getLength(); 

  bool converged = false; 
  nat ctr = 0; 
  
  auto branch2PairOptNrd2 = std::unordered_map<BranchPlain, std::array<double,3>>{}; 

  while (not converged)
    {
      converged = true; 
      for(auto b : involvedBranches)
	{
	  auto bl  = traln.getBranch(b, param); 
	  auto optTuple = GibbsProposal::optimiseBranch( traln, bl, eval, 30, param);

	  auto lengthBefore = bl.getInterpretedLength(traln, param); 
	  bl.setLength(optTuple[0]); 
	  auto lengthAfter = bl.getInterpretedLength(traln, param); 
	  auto blRatio = fabs(log(lengthBefore) - log(lengthAfter)); 

	  // tout << MAX_SCI_PRECISION << "[" << ctr << "] " << bl.toPlain() << "\t" << lengthAfter  << "\t" << SOME_FIXED_PRECISION << optTuple[1]<< MAX_SCI_PRECISION  << "\t" << optTuple[2] << std::endl; 
	  branch2PairOptNrd2[bl.toPlain()] = {{ bl.getInterpretedLength(traln, param), optTuple[1], optTuple[2] }}; 

	  converged &= blRatio < thresh || optTuple[1] > 1e-1 ;   
	  traln.setBranch(bl, param); 
	}
      ++ctr;
      // tout << "-" << std::endl;  
    }

#ifdef VERBOSE_INFO
  // TEST
  double optDiff = 0 ; 
  if(forward)
    {
      for(auto elem :  branch2PairOptNrd2) 
	{
	  auto b = elem.first.toBlDummy(); 
	  b.setConvertedInternalLength(traln, param,elem.second[0]); 
	  traln.setBranch(b,param);
	}
      eval.evaluate(traln, branch2PairOptNrd2.begin()->first, true); 
      double lnl = traln.getTr()->likelihood; 
      optDiff = lnl - initLnl; 
    }
  // END
#endif
  
  // propose in forward case, score original branches in backwards move 
  if(forward)
    {
      // propose 
      for(auto b : branch2PairOptNrd2)
	{
	  auto nrOpt = b.second[0]; 
	  auto nrd1 = b.second[1]; 
	  auto nrd2 = b.second[2]; 
	  
	  auto proposalResult = GibbsProposal::propose(nrOpt, nrd1, nrd2, rand); 

	  auto newB = b.first.toBlsDummy(); 
	  auto tmp =  b.first.toBlDummy();

	  assert( params.size() == 1); 

	  tmp.setConvertedInternalLength(traln, param, proposalResult[0]); 

	  newB.setLengths( { tmp.getLength() } );
	  
	  if(not BoundsChecker::checkBranch(newB))
	    BoundsChecker::correctBranch(newB); 
	  
	  // nasty, but necessary, because 

	  double newRealLength = 0; 
	  {
	    auto len = newB.getLength(param); 
	    auto bl = newB.toBlDummy(); 
	    bl.setLength(len); 
	    newRealLength = bl.getInterpretedLength(traln, param);
	  }

	  result.push_back(newB);
	  auto probPart = Density::lnGamma(newRealLength, proposalResult[1], proposalResult[2]); 

#ifdef VERBOSE_INFO
	  tout << MAX_SCI_PRECISION << "opt=" <<  nrOpt << "\tprop="<< newRealLength << "\t" << SOME_FIXED_PRECISION <<  log(nrOpt) -  log(newRealLength)  << "\t" << probPart << std::endl; 
#endif

	  probability *= probPart; 
	}
    }
  else 
    {
      for(auto elem : branch2PairOptNrd2)
	{
	  auto bl = elem.first.toBlDummy();
	  bl.setLength(branch2lengthOrig[elem.first]); 
	  
	  auto nrOpt = elem.second[0]; 
	  auto nrd1 = elem.second[1]; 
	  auto nrd2 = elem.second[2]; 
	  auto params = GibbsProposal::getGammaParams(nrOpt, nrd1, nrd2); 
	  auto origLen = bl.getInterpretedLength(traln,param); 

	  auto probPart = Density::lnGamma(origLen, params.first, params.second); 

	  probability *= probPart; 
#ifdef VERBOSE_INFO
	  tout << MAX_SCI_PRECISION << "opt=" << nrOpt << "\tbl="<< origLen << "\tp=" << SOME_FIXED_PRECISION << probPart << std::endl; 
#endif
	}
    }


  double ratio = 0; 
  for(auto elem : branch2PairOptNrd2)
    {
      auto bl = elem.first.toBlDummy();
      auto before = elem.second[0]; 
      bl.setLength(branch2lengthOrig[elem.first]); 
      auto after = bl.getInterpretedLength(traln, param); 
      
      double impactHere = fabs(log(after) - log(before)) ; 
      if(ratio < impactHere )
	ratio = impactHere; 
      
      // tout << MAX_SCI_PRECISION << before << "\t" << after << "\t" << impactHere << std::endl; 
    }
  
  // actually, we need to reset a bit more ... 
  for(auto elem : branch2lengthOrig)
    {
      auto b = elem.first.toBlDummy(); 
      b.setLength(elem.second); 
      traln.setBranch(b, param);
    }
  revertTree(traln,params, eval, proposeOuter );   

#ifdef VERBOSE_INFO
  if(forward)
    tout << "OPT: " << optDiff << std::endl; 
#endif

  assert(ratio >= 0); 

  return make_tuple(result, probability, ratio); 
}


SprMove SprMove::getInverseMove(const TreeAln &traln, const std::vector<AbstractParameter*> &params) const 
{
  auto result = SprMove{}; 
  auto newPath = getPathAfterMove(path);
  
  newPath.saveBranchLengthsPath(traln, params);
  
  result.setPath(newPath); 
  return result; 
}


// #define DO_INTEGRATE 

void SprMove::integrateBranches( TreeAln &traln, const std::vector<AbstractParameter*> blParams, LikelihoodEvaluator &eval, double &hastings ) const 
{
  double prevLnl = traln.getTrHandle().likelihood;  

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
      eval.evaluate(traln, b, true); 

      auto bCpy = b.toBlDummy(); 
      auto resultTmp = GibbsProposal::optimiseBranch(traln, bl, eval,  30, blParam); 
      auto result = resultTmp[0]; 

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
      eval.evaluate(traln, b, true); 
      auto resultLength = GibbsProposal::optimiseBranch(traln, bl , eval, 30, blParam)[0]; 

      bl.setLength(resultLength); 
      traln.setBranch(bl, blParam); 
      eval.evaluate(traln, b, true); 
      resultVec.push_back(bl.getInterpretedLength(traln, blParam)); 
#endif
    }

  eval.evaluate(traln, prunedSubtree, true); 
  double newLnl = traln.getTrHandle().likelihood; 

  revertTree(traln, blParams); 
  hastings = oldHastings; 

  auto bip= path.at(0).getThirdBranch(traln, path.at(1)).getBipartition(traln);
  tout << "PARS\t" << getNniDistance() << "\t" << MAX_SCI_PRECISION <<  newLnl - prevLnl << "\t"; 
  for(auto &v : resultVec)
    tout << v << "\t"; 
  tout << bip << std::endl; 
} 


std::vector<nat> SprMove::getDirtyNodes(const TreeAln &traln, bool considerOuter) const 
{
  auto result = std::vector<nat>(); 
  result.reserve(path.getNumberOfNodes()); 
  for(int i = 1; i < path.getNumberOfNodes() -1 ; ++i)
    result.push_back(path.getNthNodeInPath(i)) ;

  if(considerOuter)
    {
      auto n1 = path.getNthNodeInPath(0); 
      auto n2 = path.getNthNodeInPath(path.getNumberOfNodes()-1); 
      if(not traln.isTipNode(n1))
	result.push_back(n1); 
      if(not traln.isTipNode(n2))
	result.push_back(n2); 
    }
  
  return result; 
} 


bool SprMove::returnCommonBranchesAfter(const BranchPlain &branch,  BranchPlain &result) const
{
  // not on path or pruned subtree  
  if (not ( path.nodeIsOnPath(branch.getPrimNode() ) && path.nodeIsOnPath(branch.getSecNode() )) ) 
    {
      result = branch; 
      return true; 
    }

  // first branch  
  auto firstNode = path.getNthNodeInPath(0); 
  if( branch.hasNode(firstNode) )
    {
      result = BranchPlain(firstNode, path.getNthNodeInPath(2)); 
      return true; 
    }
  
  // very last branch 
  auto lastnode = path.getNthNodeInPath(path.getNumberOfNodes() - 1 ); 
  if(branch.hasNode(lastnode))
    {
      result = BranchPlain(lastnode, path.getNthNodeInPath(1)); 
      return true; 
    }

  return false; 
}


BranchPlain SprMove::mapBranchNniStepsAfter(const BranchPlain &branch)const  
{
  auto result = BranchPlain(); 
  if(returnCommonBranchesAfter(branch, result))
    return result; 
  
  // second to last needs special attention 
  auto secToLast = path.getNthNodeInPath(path.getNumberOfNodes() - 2 ); 
  if(branch.hasNode(secToLast))
    return BranchPlain(secToLast, path.getNthNodeInPath(1)); 
  
  // else, map shift the branches one towards the insertion position
  for(nat i = 0; i < path.size() ; ++i)
    {
      if(branch.equalsUndirected(path.at(i))) 
	return path.at(i+1); 
    }

  // we must have found something 
  assert(0); 
  return BranchPlain(); 
}



BranchPlain SprMove::mapBranchSingleMapAfter(const BranchPlain &branch ) const
{
  auto result = BranchPlain(); 
  if(returnCommonBranchesAfter(branch, result))
    return result; 
  
  // map the mapped branch  
  if(path.at(1).equalsUndirected(branch)) 
    return BranchPlain(path.getNthNodeInPath(1), path.getNthNodeInPath(path.getNumberOfNodes()-2)); 
  
  // otherwise, everything stays the same 
  return branch; 
} 


template<typename T1,typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1,T2> &elem)
{
  out << elem.first << "," << elem.second ; 
  return out; 
}

static void nameSubtree(const TreeAln &traln, const BranchPlain &b, std::string name,  branch2PairNameNum &result, nat depth )
{
  auto elem = std::make_pair(name,depth); 
  assert(result.find(b) == result.end()); 
  result[b] = elem; 

  if(not b.isTipBranch(traln))
    {
      auto desc = traln.getDescendents(b.getInverted()); 
      nameSubtree(traln, desc.first, name, result, depth+ 1);
      nameSubtree(traln, desc.second, name, result, depth+ 1);
    }
}

// must be applied before the move!  
branch2PairNameNum SprMove::getNames(const TreeAln &traln, bool isNni ) const 
{
  auto result = branch2PairNameNum{}; 
  auto subTreeNames = std::vector<std::string> { "A","B", "C", "D" ,"E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q" , "R", "S", "T", "U", "V", "W", "X", "Y", "Z"}; 
  
  nat trCtr = 0;
  
  auto start = BranchPlain(path.getNthNodeInPath(1), path.getNthNodeInPath(0)); 
  nameSubtree(traln,start , subTreeNames[trCtr], result, 1); 
  ++trCtr; 

  for(nat i = 1; i < path.size(); ++i)
    {
      auto start = path.at(i-1).getThirdBranch(traln,path.at(i)); 
      nameSubtree(traln, 
		  start, 
		  subTreeNames[trCtr], result, 1); 
      ++trCtr; 
      
      auto branch = path.at(i); 
      if(i == path.size() - 1 )	//  last branch (is a new subtree again )
	{
	  nat numNodes = path.getNumberOfNodes(); 
	  nameSubtree(traln, BranchPlain(path.getNthNodeInPath(numNodes-2), path.getNthNodeInPath(numNodes-1)), subTreeNames[trCtr], result, 1); 
	  ++trCtr; 	    
	}
      else if(isNni)
	{
	  assert(result.find(branch) == result.end()) ; 
	  auto elem = make_pair("INNER", i-1); 
	  result[branch] = elem ; 
	}
      else if(not isNni)
	{
	  assert(result.find(branch) == result.end()) ; 
	  auto elem = make_pair("a", 0); 
	  if( i == 1 )
	    elem = make_pair("SWITCH", 0); 
	  else 
	    elem = make_pair("INNER", i-2); 
	  result[branch] = elem; 
	}
      else 
	assert(0); 
    }

  assert(result.size() ==  traln.getNumberOfBranches() ) ; 
  return result; 
}



// for debug 
void SprMove::printBothMappings() const 
{
  tout << "for path (NNI,Single):\t"  << path  << std::endl; 
  for(nat i = 0; i < path.size() ; ++i)
    {
      auto branch = path.at(i); 
      tout << branch << "\t" << mapBranchNniStepsAfter(branch) << "\t" << mapBranchSingleMapAfter(branch) << std::endl; 
    }
}


std::vector<SprMove> SprMove::getAllUniqueMoves(const TreeAln& traln, nat dist)
{
  auto result = std::vector<SprMove>{}; 
  
  assert(dist > 0); 

  for(auto pruneBranch : traln.extractBranches())  
    {
      auto pP = pruneBranch.getPrimNode();  

      if(not traln.isTipNode( pruneBranch.getPrimNode()))
	{
	  auto distBranches = traln.getBranchesByDistance(pruneBranch, dist + 1, false ); 

	  for(auto insertBranch : distBranches)
	    {
	      auto iP = insertBranch.getPrimNode(),
		iS = insertBranch.getSecNode(); 

	      if(dist == 1 
		 &&   
		 ( ( BranchPlain(pP, iP ).exists(traln) && pP > iP  ) 
		   || ( BranchPlain(pP, iS ).exists(traln) && pP > iS  )) )
		{
		  // tout << "skipping " << pruneBranch << "," << insertBranch << std::endl; 
		  continue; 
		}

	      auto elem = SprMove{}; 
	      elem.extractBranchesOnly(traln, pruneBranch, insertBranch, elem.getPathHandle()); 
	      result.push_back( elem);
	    }
	}
      
      auto inversePruneBranch = pruneBranch.getInverted(); 
      pP = inversePruneBranch.getPrimNode(); 

      if(not traln.isTipNode(inversePruneBranch.getPrimNode()))
	{
	  auto distBranches = traln.getBranchesByDistance(inversePruneBranch, dist + 1, false); 
	  
	  for(auto insertBranch: distBranches)
	    {
	      auto iP =  insertBranch.getPrimNode(),
		iS = insertBranch.getSecNode(); 

	      if(dist == 1
		 && (  ( BranchPlain(pP, iP).exists(traln) && pP > iP  ) 
		       || (BranchPlain(pP, iS).exists(traln) && pP > iS ))) 
		{
		  // tout << "skipping " << pruneBranch << "," << insertBranch << std::endl; 
		  continue; 
		}

	      auto elem = SprMove{}; 
	      elem.extractBranchesOnly(traln, inversePruneBranch, insertBranch, elem.getPathHandle()); 
	      result.push_back(elem );
	    }
	}
    }
  
  return result; 
}


