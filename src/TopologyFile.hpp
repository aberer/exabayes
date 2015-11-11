#ifndef TOPOLOGY_FILE
#define TOPOLOGY_FILE

#include <sstream>
#include <string>
#include <iostream> 

class TopologyFile
{
public: 
  TopologyFile(std::string workdir, std::string runname, nat runid, nat couplingId)
  { 
    std::stringstream ss ; 
    // TODO portability 
    ss << workdir <<  ( workdir.compare("") == 0  ? "" :  "/" )  << "ExaBayes_topologies." << runname << "." << runid ; 
    if(couplingId != 0 )
      ss << ".hot-"<<  couplingId; 
    fullFilename = ss.str();

    std::ofstream fh (fullFilename); 
    fh << "" ; 
    fh.close(); 
  }


  void initialize(const TreeAln& traln, nat someId) const 
  {
    std::ofstream fh(fullFilename,std::fstream::app|std::fstream::out ); 
    fh << "#NEXUS" << std::endl
       << "[ID: " << someId <<  "]" << std::endl
       << "[Param: tree]" << std::endl 
       <<  "begin trees;" << std::endl
       << "\ttranslate" <<  std::endl; 
    
    bool isList = false; 
    nat numTax = traln.getNumberOfTaxa(); 
    for(nat i = 0; i < numTax; ++i)
      {
	isList = ( i == numTax - 1 )  ; 
	fh << "\t" <<  i+1 << "\t"
	      << traln.getTr()->nameList[i+1]  << (isList ? ";" : ",")<< std::endl; 
      }    
    fh.close(); 
  }  

  
  void sample(const TreeAln &traln, nat gen) const 
  {    
    std::ofstream fh(fullFilename,std::fstream::app|std::fstream::out ); 
    TreePrinter tp(true, false, false);
    std::string treeString = tp.printTree(traln);
    fh << "\ttree gen." << gen << " = [&U] " << treeString << std::endl; 
    fh.close(); 
  }

  
  void finalize() const 
  {
    std::ofstream fh(fullFilename, std::fstream::app|std::fstream::out ); 
    fh << "end;" << std::endl; 
    fh.close(); 
  }

private: 
  std::string fullFilename; 
}; 

#endif
