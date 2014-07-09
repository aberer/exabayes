#ifndef _ABSTRACTPROPOSAL_H
#define _ABSTRACTPROPOSAL_H

#include <string>
#include <vector>

class OptimizedParameter;  

#include "model/Category.hpp"
#include "mcmc/SuccessCounter.hpp" 
#include "math/Randomness.hpp"
#include "priors/PriorBelief.hpp"
#include "system/GlobalVariables.hpp"
#include "eval/LikelihoodEvaluator.hpp"
#include "TreeRandomizer.hpp"
#include "system/Serializable.hpp"
#include "DistributionProposer.hpp"
#include "GammaProposer.hpp"



class AbstractProposal : public Serializable
{
public: 
  AbstractProposal( Category cat, std::string  _name, double weight, double minTuning, double maxTuning, bool needsFullTraversal )  ; 
  AbstractProposal( const AbstractProposal& rhs)  ;   
  virtual ~AbstractProposal(){}
  AbstractProposal& operator=(const AbstractProposal &rhs) = delete;  

  bool isUsingOptimizedBranches() const {return _usingOptimizedBranches; } 

  std::array<bool,3> getBranchProposalMode() const ; 

  // you MUST implement all virtual methods in your derived
  // proposal. Here, the signatures are set to 0, this must not
  // be the case in the derived proposal. This 0 keyword makes it
  // impossible to create an instance of AbstractProposal and forces
  // you to implement these methods, when you derive from
  // AbstractProposal.

  /**
     @brief determines the proposal, applies it to the tree / model, updates prior and hastings ratio
   */ 
  virtual void applyToState(TreeAln &traln, PriorBelief &prior, log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval) = 0; 
  /**
     @brief evaluates the proposal 
     @todo remove the prior, we should not need it here 
   */ 
  virtual void evaluateProposal(LikelihoodEvaluator &evaluator, TreeAln &traln, const BranchPlain &branchSuggestion) = 0; 
  /** 
      @brief resets the tree to its previous state; corrects the prior, if necessary (@todo is this the case?)
   */ 
  virtual void resetState(TreeAln &traln) = 0 ; 
  /** 
      @brief tunes proposal parameters, if available  
   */ 
  virtual void autotune() = 0  ;
  
  virtual AbstractProposal* clone() const = 0;  
  /** 
      @brief gets the relative weight of this proposal 
   */ 
  double getRelativeWeight() const { return _relativeWeight; }
  /**
     @brief sets the relative weight of this proposal 
   */ 
  void setRelativeWeight(double tmp) { _relativeWeight = tmp; }
  
  virtual BranchPlain determinePrimeBranch(const TreeAln &traln, Randomness& rand) const = 0; 
  /** 
      @brief gets the category 
   */ 
  Category getCategory() const {return _category; }
  /** 
      @brief gets the name of the proposal 
   */ 
  std::string getName() const {return _name; }
  /** 
      @brief gets nodes that are invalid by executed the proposal 
   */ 
  virtual std::vector<nat> getInvalidatedNodes(const TreeAln &traln) const = 0;  
  /** 
      @brief inform proposal about acceptance
   */ 
  virtual void accept() {_sctr.accept();}
  /**
     @brief inform proposal about rejection 
  */ 
  virtual void reject() {_sctr.reject();}
  /** 
      @brief gets the success counter 
   */ 
  const SuccessCounter& getSCtr()  const { return _sctr; }
  /** 
      @brief gets the number of proposal invocation, since it has been tune the last time 
   */ 
  int getNumCallSinceTuning() const { return _sctr.getRecentlySeen(); }
  /** 
      @brief add a parameter to be integrated over to the proposal 
   */ 
  void addPrimaryParameter(std::unique_ptr<AbstractParameter> var) {_primaryParameters.push_back(std::move(var)) ; }
  /** 
      @brief add a parameter  that is integrated over as a by-product of this proposal 
   */ 
  void addSecondaryParameter(std::unique_ptr<AbstractParameter> var) {_secondaryParameters.push_back(std::move(var)) ; }
  /** 
      @brief indicates whether this proposal needs a full traversal 
   */ 
  bool isNeedsFullTraversal() const {return _needsFullTraversal; }
  std::vector<AbstractParameter*> getBranchLengthsParameterView() const ; 

  /** 
      @brief update the hatsings
      @param valToAdd is the absolute proposal ratio that shall be added 
   */ 
  // static void updateHastingsLog(double &hastings, double valToAdd, std::string whoDoneIt); 

  std::vector<nat> getAffectedPartitions() const ; 

  std::ostream& printNamePartitions(std::ostream &out); 
  std::ostream& printShort(std::ostream &out)  const ;

  std::vector<AbstractParameter*> getPrimaryParameterView() const; 
  std::vector<AbstractParameter*> getSecondaryParameterView() const ; 
  /** 
      @brief prepare for set execution (only relevant for branch length + node slider)
   */ 
  virtual std::pair<BranchPlain,BranchPlain> prepareForSetExecution(TreeAln &traln, Randomness &rand)  = 0; 

  // we need to implement these 
  virtual void serialize( std::ostream &out)  const;  
  virtual void deserialize( std::istream &in ) ; 
  
  void setPreparedBranch(BranchPlain b ) {_preparedBranch = b;  }
  void setOtherPreparedBranch(BranchPlain b){_preparedOtherBranch = b; }

  void setInSetExecution(bool exec) { _inSetExecution = exec;  }
  void setId(nat id){_id = id; }
  nat getId() const {return _id; }

  /** 
      @brief writes proposal specific (tuned) parameters
   */ 
  virtual void writeToCheckpointCore(std::ostream &out)const  = 0 ;  
  /** 
      @brief reads proposal specific (tuned) parameters 
   */ 
  virtual void readFromCheckpointCore(std::istream &in) = 0; 

  virtual void prepareForSetEvaluation( TreeAln &traln, LikelihoodEvaluator& eval) const  {} 


  virtual void printParams(std::ostream &out)  const {} 

  friend std::ostream&  operator<< ( std::ostream& out , const AbstractProposal& rhs); 

  void setForProteinsOnly(bool val)  {_forProteinOnly = val; }
  bool isForProteinOnly() const {return _forProteinOnly; }

  double tuneParameter(int batch, double accRatio, double parameter, bool inverse ); 

  nat getNumTaxNeeded() const {return _numTaxNeeded; }

  virtual void extractProposer(TreeAln& traln, const OptimizedParameter& param) { }

protected:   
  std::string _name;   
  SuccessCounter _sctr; 
  Category _category; 
  std::vector<std::unique_ptr<AbstractParameter> > _primaryParameters; // it is the  primary purpose of this proposal to integrate over these parameters (in most cases only 1) 
  std::vector<std::unique_ptr<AbstractParameter> > _secondaryParameters;  // as a by-product also these random variables are changed 
  double _relativeWeight; 
  bool _needsFullTraversal; 
  bool _inSetExecution;

  // meh 
  BranchPlain _preparedBranch; 
  BranchPlain _preparedOtherBranch; 
  
  nat _id; 

  bool _forProteinOnly; 

  double _minTuning; 
  double _maxTuning; 

  nat _numTaxNeeded; 
  bool _usingOptimizedBranches; 
}; 

#endif

