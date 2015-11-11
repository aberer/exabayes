/** 
    @file SampleMaster.hpp
    
    @brief represents various chains sampling the posterior probability space
    
    Despite of its modest name, this is in fact the master class.  
 */ 

#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <vector>

#include "proposals/ProposalSet.hpp"
#include "config/BlockRunParameters.hpp"
#include "config/BlockProposalConfig.hpp"
#include "config/CommandLine.hpp"
#include "CoupledChains.hpp"
#include "config/ConfigReader.hpp"
#include "comm/ParallelSetup.hpp"
#include "system/time.hpp"
#include "system/Serializable.hpp"
#include "file/DiagnosticsFile.hpp"


class SampleMaster : public Serializable
{
public:   
  SampleMaster() ; 

  /** 
      @brief initializes the runs  
      @notice this is the top-level function 
   */ 
  void initializeRuns(Randomness rand); 
  // nat peekNumTax(std::string filePath); 
  /** 
      @brief cleanup, once finished
   */ 
  void finalizeRuns();  
  /** 
      @brief starts the MCMC sampling   
   */ 
  void run(); 
  /** 
      @brief initializes the config file 
   */ 
  std::tuple<std::vector<std::unique_ptr<AbstractParameter> > , std::vector<std::unique_ptr<AbstractProposal> > , std::vector<ProposalSet> >  
  processConfigFile(string configFileName, const TreeAln &tralnPtr ) ; 
  void initializeWithParamInitValues(TreeAln &tree , const std::vector<AbstractParameter*> &params , bool hasBl ) const ; 
  /** 
      @brief EXPERIMENTAL 
   */ 
  void branchLengthsIntegration(Randomness &rand)  ;  
  /** 
      @brief print information about the alignment  
   */ 
  void printAlignmentInfo(const TreeAln &traln);   
  /** 
      @brief gets the runs 
   */ 
  const std::vector<CoupledChains>& getRuns() const {return _runs; }
  
  virtual void deserialize( std::istream &in ) ; 
  virtual void serialize( std::ostream &out) const ;

  void setCommandLine(CommandLine cl) { _cl = cl; }

  void setParallelSetup(ParallelSetup* pl ) { _plPtr = pl; }

  void synchronize( CommFlag flags )  ; 

private: 
  void printParameters(const TreeAln &traln, const std::vector<unique_ptr<AbstractParameter> > &params) const; 
  void printProposals(const std::vector<unique_ptr<AbstractProposal> > &proposals, const std::vector<ProposalSet> &proposalSets  ) const ; 
  void printInitializedFiles() const; 
  std::vector<std::string> getStartingTreeStrings(); 
  void informPrint(); 
  void printInitialState(); 
  LikelihoodEvaluator createEvaluatorPrototype(const TreeAln &initTree, bool useSEV); 
  void writeCheckpointMaster(); 
  void initializeFromCheckpoint(); 
  std::pair<double,double> convergenceDiagnostic(nat &begin, nat &end); 
  std::vector<bool> initTrees(std::vector<TreeAln> &trees, randCtr_t seed, std::vector<std::string> startingTreeStrings, const std::vector<AbstractParameter*> &params); 
  bool initializeTree(TreeAln &traln, std::string startingTree, Randomness &treeRandomness, const std::vector<AbstractParameter*> &params); 
  CLOCK::system_clock::time_point printDuringRun(nat gen) ; 
  std::string getOrCreateBinaryFile() const ; 
  void catchRunErrors() const ; 

private:			// ATTRIBUTES 
  std::vector<CoupledChains> _runs; 
  ParallelSetup* _plPtr;  	// non-owning!
  CLOCK::system_clock::time_point _initTime; 
  BlockRunParameters _runParams;  
  CommandLine _cl; 
  CLOCK::system_clock::time_point _lastPrintTime; 
  DiagnosticsFile _diagFile; 
};  

#endif
