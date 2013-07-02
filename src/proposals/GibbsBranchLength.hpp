#include "GibbsProposal.hpp" 
#include "BranchLengthMultiplier.hpp"

class GibbsBranchLength : public BranchLengthMultiplier
{
public: 
  GibbsBranchLength(LikelihoodEvaluatorPtr _eval)
    : BranchLengthMultiplier(0)
    , eval(_eval)
  {
    name = "estGibbsBL"; 
    category = Category::BRANCH_LENGTHS; 
    relativeWeight = 20 ; 
  } 

  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
  {
    Branch  b = proposeBranch(traln, rand);     
    b.setLength(b.findNodePtr(traln)->z[0]); 
    double initBl = b.getLength(); 


    nodeptr p = b.findNodePtr(traln); 
    savedBranch = b;     

    
    GibbsProposal::drawFromEsitmatedPosterior(b, eval, traln, rand, initBl, 5, hastings); 
    double newZ = b.getLength();
    traln.setBranchLengthBounded(newZ, 0, b.findNodePtr(traln)); 

    auto brPr = primVar[0]->getPrior();
    prior.updateBranchLengthPrior(traln,initBl, newZ, brPr);
  }

  virtual ~GibbsBranchLength(){}  
  virtual void autotune(){}

  virtual AbstractProposal* clone() const  { return new GibbsBranchLength(*this); }  
  
private: 
  LikelihoodEvaluatorPtr eval; 

}; 
