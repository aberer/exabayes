#ifndef _MYTEMPLATEPROPOSAL_H  	// use whatever you want as an include guard (must NOT occur twice! )
#define _MYTEMPLATEPROPOSAL_H

#include "AbstractProposal.hpp"


// if you created a new source file (=> .cpp), then inform automake
// about it. that's done with the
// make my-update-src
// command. Probably it will run configure again. 
// If for any reason nothing works any more after it, have a look, if the source files are correct in 
// src.am (you can also manually modify this source, s.t. it works). 


// if you do not use the namespace below, you have to
// add the std:: prefix to many things (e.g., the vector is std::vector then). 
// using the qualifier is better practice, but inconvenient. 
// using namespace std; 


class MyTemplateProposal : public AbstractProposal
{
public: 
  MyTemplateProposal( double aVariable);

  // everything copied over from AbstractProposal that has a = 0
  // there. This declares the methods pure virtual, and the derived
  // proposal HAS to implement theese methods.
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand, LikelihoodEvaluator& eval) ; 
  virtual void evaluateProposal(  LikelihoodEvaluator &evaluator, TreeAln &traln) ; 
  virtual void resetState(TreeAln &traln); 
  virtual void autotune()  ;
  virtual AbstractProposal* clone() const ;  

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return Branch(0,0);}
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}
  
  // if your proposal has  parameters that are tuned   
  // this is very straight forward, just check out the other proposals 
  virtual void readFromCheckpointCore(std::istream &in); 
  virtual void writeToCheckpointCore(std::ostream &out) const; 


private: 
  void privateMethod(TreeAln &traln, Randomness& rand)  ; 	// some helper method. It is in the private sector, s.t. it cannot be accessed from the outside 

  double aVariable; 		// needed for remembrance  
  double aTuningParamer; 	// for instance a multiplier we draw from  

};


#endif
