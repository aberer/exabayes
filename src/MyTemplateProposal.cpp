#include "MyTemplateProposal.hpp"



// we need a constructor!  
MyTemplateProposal::MyTemplateProposal(double _relativeProbability, double aVariable)
{
  this->relativeProbability = _relativeProbability; //  the constructor argument should not have the exact same name as the member variable

  name = "MyTemplateProposal"; 	
  category = TOPOLOGY ; 	// check out categoryType.h
  
  // ptype = E_TBR; 		// actually not used any more... 

}




void MyTemplateProposal::applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand)
{
  tree *tr = traln.getTr(); 	// our beloved tree
  
  // just writing nonsense here, s.t. the compiler does not complain 
  tr->likelihood = 0; 

  // do some fancy stuff that basically is the proposal. 

  // NOTICE: the argument is double &hastings. This is a specific c++
  // construct called reference. It behaves like a pointer (i.e., if
  // you modify this value, the hastings value OUTSIDE the function
  // will also be modified), but syntactically it is handled like a
  // raw value (no de-referencing needed) => hastings = 1; works!

  // for changing the hastings, use this function: 
  double somethingToAddToTheHastings = 3 ; 
  updateHastings(hastings, somethingToAddToTheHastings , name);
  // internally the hastings is on log-scale. Your value to be added
  // should NOT be on the log-scale

  
  // since rand is an argument, we easily can draw random numbers:
  double r =  rand.drawRandDouble01();
  cout << "the value was "<<  r << endl; 


  // maybe we need our private method as a helper
  privateMethod(traln,rand);


  // assume we have changed the branch lengths. Then we have to modify the prior as well.   
  double oldBL = 0; 		// the REAL vaule 
  double newBL = 0; 	
  prior.updateBranchLength(oldBL, newBL); // updates the prior, also
					  // check out the other
					  // methods for modification.

  // for verifying that your prior modifications are correct, check
  // out the DEBUG_LNPR in common.h
}


void MyTemplateProposal::evaluateProposal(TreeAln &traln, PriorBelief &prior) 
{
  // somehow determine where what should be evaluated
  // the prior is NOT changed here
}



void MyTemplateProposal::resetState(TreeAln &traln, PriorBelief &prior)
{  
  // reset your modifications to the tree   
  
  // dont forget to update the prior: 
  double newBL = 0; 		// we probably have stored that as a private member variable
  double oldBL = 0; 
  
  prior.updateBranchLength(newBL, oldBL); // values are switched this time  
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
  return new MyTemplateProposal(relativeProbability, aVariable);
}








void MyTemplateProposal::privateMethod(TreeAln& traln, Randomness &rand)
{
  cout << "i am modifying the tree is some manner" << endl; 
}