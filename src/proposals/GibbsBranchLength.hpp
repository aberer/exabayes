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


    assert(0);
    // TODO 
    // b.setLength(b.findNodePtr(traln)->z[0]); 
    // double initBl = b.getLength( ); 
    double initBl = 0; 
    
    assert(primaryParameters.size() == 1); 
    savedBranch = b; 
    
    auto param = primaryParameters[0].get();
    GibbsProposal::drawFromEsitmatedPosterior(b, *eval, traln, rand,  MAX_ITER, hastings, param); 
    double newZ = b.getLength(param);
    traln.setBranch(b, getSecondaryParameterView());    

    prior.updateBranchLengthPrior(traln,initBl, newZ, param);  
  }

  virtual void autotune(){}

  // virtual Branch prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return Branch(0,0);}
  virtual std::pair<Branch,Branch> prepareForSetExecution(TreeAln &traln, Randomness &rand)  { return std::pair<Branch, Branch> (Branch(0,0),Branch(0,0) );}

  virtual AbstractProposal* clone() const  { return new GibbsBranchLength(*this); }  

private: 
  std::unique_ptr<LikelihoodEvaluator> eval; 
}; 
