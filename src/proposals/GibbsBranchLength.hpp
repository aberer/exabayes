#include "GibbsProposal.hpp" 
#include "BranchLengthMultiplier.hpp"
#include "Category.hpp"


#define MAX_ITER 30 

class GibbsBranchLength : public BranchLengthMultiplier
{
public: 
  GibbsBranchLength(std::unique_ptr<LikelihoodEvaluator> _eval)
    : BranchLengthMultiplier(0)
    , eval(std::move(_eval))
  {
    name = "estGibbsBL"; 
    category = Category::BRANCH_LENGTHS; 
    relativeWeight = 20 ; 
  } 

  GibbsBranchLength(const GibbsBranchLength& rhs)
    : BranchLengthMultiplier(rhs)
    , eval(rhs.eval->clone())
  {
  }



  virtual void applyToState(TreeAln &traln, PriorBelief &prior, double &hastings, Randomness &rand) 
  {
    Branch b = proposeBranch(traln, rand);     
    b.setLength(b.findNodePtr(traln)->z[0]); 
    double initBl = b.getLength(); 

    savedBranch = b;     
    
    GibbsProposal::drawFromEsitmatedPosterior(b, *eval, traln, rand,  MAX_ITER, hastings); 
    double newZ = b.getLength();
    traln.setBranch(b);    

    auto brPr = primaryParameters[0]->getPrior();
    prior.updateBranchLengthPrior(traln,initBl, newZ, brPr); // 
  }

  virtual void autotune(){}

  virtual AbstractProposal* clone() const  { return new GibbsBranchLength(*this); }  

private: 
  std::unique_ptr<LikelihoodEvaluator> eval; 
}; 
