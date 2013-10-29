#include "TreeRandomizer.hpp"

#include <vector>
#include "Branch.hpp"

void TreeRandomizer::createParsimonyTree(TreeAln &traln, Randomness& rand)
{
  nat r = rand();  

  traln.unlinkTree();
  traln.getTr()->randomNumberSeed = r; 

#if HAVE_PLL != 0
  makeParsimonyTreeFast(traln.getTr(), traln.getPartitionsPtr());
#else 
  makeParsimonyTreeFast(traln.getTr()); 
#endif
}


void TreeRandomizer::randomizeTree(TreeAln &traln, Randomness& rand )
{
  for(nat i = 1 ; i < traln.getNumberOfNodes() + 1 ; ++i)
    {
      auto p = traln.getNode(i);
      p->back = NULL; 
      p->next->back = NULL; 
      p->next->next->back = NULL;       
    }

  // start with the simple tree 
  auto a = traln.getNode(1 ),
    b = traln.getNode(2 ), 
    c = traln.getNode( 3 ), 
    inner = traln.getNode(traln.getNumberOfTaxa() + 1); 
  
  traln.clipNodeDefault(inner,a); 
  traln.clipNodeDefault(inner->next, b); 
  traln.clipNodeDefault(inner->next->next,c); 

  for(nat i = 4; i < traln.getNumberOfTaxa() +1 ;  ++i)
    {
      inner = traln.getNode(traln.getNumberOfTaxa() + i-2 );       
      auto taxonP = traln.getNode(i); 
      traln.clipNodeDefault(taxonP, inner); 
      
      auto b = drawBranchUniform_helper(traln, rand, i-1);

      auto p1 = b.findNodePtr(traln),
	p2 = p1->back; 
      
      traln.clipNodeDefault(p1, inner->next); 
      traln.clipNodeDefault(p2, inner->next->next); 
    }
 
}


BranchPlain TreeRandomizer::drawInnerBranchUniform( const TreeAln& traln, Randomness &rand)  
{
  bool acc = false;   
  int node = 0; 
  nodeptr p = nullptr; 
  while(not acc)
    {      
      node = drawInnerNode(traln, rand); 
      p = traln.getNode(node); 
      
      nat numTips = 0; 
      if(  traln.isTipNode(p->back) ) 
	numTips++; 
      if(traln.isTipNode(p->next->back))
	numTips++;
      if(traln.isTipNode(p->next->next->back))
	numTips++; 
      
      assert(numTips != 3); 
      
      acc = numTips == 0 || rand.drawRandDouble01() <  (3. - double(numTips)) / 3.;       
    }
  assert(node != 0); 
  
  std::vector<nat> options; 
  if(not traln.isTipNode(p->back))
    options.push_back(p->back->number); 
  if(not traln.isTipNode(p->next->back))
    options.push_back(p->next->back->number); 
  if(not traln.isTipNode(p->next->next->back))
    options.push_back(p->next->next->back->number); 

  nat other = 0;
  if(options.size() == 1 )
    other = options[0]; 
  else 
    other = options.at(rand.drawIntegerOpen(options.size()));   
  return BranchPlain(node, other); 
}


nat TreeRandomizer::drawInnerNode(const TreeAln& traln, Randomness &rand )  
{    
  nat curNumTax = traln.getNumberOfTaxa(); 
  nat res =  1 + curNumTax + rand.drawIntegerClosed(curNumTax - 3 );   
  return res; 
}


/** 
    @brief draw a branch that has an inner node as primary node   
    
    => equals draw subtree uniformly
 */ 
BranchPlain TreeRandomizer::drawBranchWithInnerNode(const TreeAln& traln,Randomness &rand)  
{
  nat idA = drawInnerNode(traln, rand); 
  nat r = rand.drawIntegerClosed(2);  
  nodeptr p = traln.getNode(idA); 
  assert(not traln.isTipNode(p)) ; 

  BranchPlain b; 
  switch(r)
    {
    case 0: 
      b = BranchPlain(idA, p->back->number); 
      break; 
    case 1: 
      b = BranchPlain(idA, p->next->back->number); 
      break; 
    case 2: 
      b = BranchPlain(idA, p->next->next->back->number); 
      break; 
    default: assert(0); 
    }
 
  return b; 
}



BranchPlain TreeRandomizer::drawBranchUniform(const TreeAln & traln, Randomness &rand )  
{
  return drawBranchUniform_helper(traln, rand, traln.getNumberOfTaxa()); 
}

BranchPlain TreeRandomizer::drawBranchUniform_helper(const TreeAln &traln, Randomness &rand , nat curNumTax)  
{ 
  int randId = 0; 
  double r = rand.drawRandDouble01(); 

  // for the randomization part, i assume that when a tree is built
  // successively, I assume that the inner nodes start with an offset
  // in the nodeptr array that is the number of trees in the final
  // tree
  if( r <= 0.75) 		// draw an inner node 
    randId = 1 + traln.getTr()->mxtips + rand.drawIntegerClosed(curNumTax -3 )  ; 
  else 				// draw a tip 
    randId = rand.drawIntegerOpen(curNumTax) + 1 ;           

  auto p = traln.getNode(randId); 
  nat thisNode = randId,
    thatNode = 0;   
  if(traln.isTipNode(p))
    {
      thatNode = p->back->number; 
    }
  else 
    {
      int r = rand.drawIntegerOpen(3); 
          switch(r)
  	{
  	case 0 : 
  	  thatNode = p->back->number; 
  	  break; 
  	case 1 : 
  	  thatNode = p->next->back->number; 
  	  break; 
  	case 2: 
  	  thatNode = p->next->next->back->number; 
  	  break; 
  	default: assert(0); 
  	}
    }

  return BranchPlain(thisNode, thatNode); 
}
