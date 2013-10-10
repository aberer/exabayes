#include "MyTemplateProposal.hpp"



// this MUST be initialized like that here. It defines the relative
// weight for all proposals of this type. This variable can be
// overridden in the config file.

// the variable MUST be defined in the header file first

// we need a constructor!  
MyTemplateProposal::MyTemplateProposal( double aVariable)
  : AbstractProposal( Category::TOPOLOGY, "MyTemplateProposal")
{

  relativeWeight = 2.0; 
}


void MyTemplateProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
  // i do not have an idea what to write here, maybe check out
  // ParameterProposal, I documented the apply-funciton there .
  
}


void MyTemplateProposal::evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) 
{
  // the previous evluation scheme has been replaced with the
  // LikelihoodEvaluator class. 

  // let's assume, we have 4 partitions and this proposal proposes for
  // 2 frequency parameters that have 2 partitions each => the scheme is (1,2/3,4)

  // if we want to evaluate all four partitions (but not necessarily
  // the entire aligment), then we'd have to do this:
  
  std::vector<nat> allRelevantPartitions; 
  for(nat i = 0; i < primaryParameters.size(); ++i)  
    {
      std::vector<nat> partitionsHere = primaryParameters[i]->getPartitions(); 
      for(nat j = 0; j < partitionsHere.size(); ++j)
	allRelevantPartitions.push_back(partitionsHere[j]); 
    }

  
  // evaluator->evaluatePartitions(traln, allRelevantPartitions, true); 
  
  
  // let's do it again and use some c++ conveniences 
  std::vector<nat> allRelevantPartitions2; 
  // range-iteration over the vector primvar. Due to the auto-keyword,
  // you do not have to remember the exact type (could be
  // difficult). However, here, you MUST add a "&". This means that
  // you do not want a copy of the parameter object, but the actual
  // object
  for(auto &var : primaryParameters)	
    for(auto p : var->getPartitions()) // iterate over all partitions in this parameter 
      allRelevantPartitions2.push_back(p); 
  // evaluator.evaluatePartitions(traln, allRelevantPartitions2,true); 

  // when in doubt, there is always 
  tree *tr = traln.getTr(); 
  auto toEval = BranchPlain(tr->start->number, tr->start->back->number); // a bit clunky 
  evaluator.evaluate(traln, toEval, true); 
}



void MyTemplateProposal::resetState(TreeAln &traln)
{  

}



void MyTemplateProposal::autotune() 
{
  // do nothing
}


// returns the base class type. This way, we can store the proposal in
// a vector<AbstractProposal*> and do not even have to know which
// proposal it is. We simply can use the methods that are implemented
// in every proposal.

// But: we need this tiny clone method. It is almost always the same
// everywhere and just constructs a new proposal (since every chain
// wants to have its own instance of a proposal)
AbstractProposal* MyTemplateProposal::clone() const 
{
  return new MyTemplateProposal( *this);
}


void MyTemplateProposal::privateMethod(TreeAln& traln, Randomness &rand)
{
  std::cout << "i am modifying the tree is some manner" << std::endl; 
}



// checkpointing stuff 
void MyTemplateProposal::readFromCheckpointCore(std::istream &in)
{  
  // assume our aTuningParamer needs to be tuned. Then, we have to
  // read it from the previous checkpointing file: 
  
  aTuningParamer = cRead<double>(in); 
}


void MyTemplateProposal::writeToCheckpointCore(std::ostream &out) const
{ 
  // and we have to write it to the filen, once we checkpoint 
  cWrite<double>(out, aTuningParamer); 	  
}
