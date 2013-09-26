#include "TopologyFile.hpp"
#include "GlobalVariables.hpp"
#include "ParallelSetup.hpp"
#include "parameters/AbstractParameter.hpp"
#include "TreePrinter.hpp"
#include "parameters/BranchLengthsParameter.hpp"

#include <cassert>

void genericExit(int code); 

TopologyFile::TopologyFile(std::string workdir, std::string runname, nat runid, nat couplingId, nat _paramNum, bool _hasManyTopoloFiles) 
  : runid(runid)
  , couplingId(couplingId)
  , paramNum(_paramNum)
  , hasManyTopoloFiles(_hasManyTopoloFiles)
{ 
  auto&& ss = std::stringstream{}; 
  
  ss <<  OutputFile::getFileBaseName(workdir) << "_topologies." << runname << "." << runid ; 
  if(hasManyTopoloFiles)
    ss << ".tree." << paramNum; 
  
  if(couplingId != 0 )
    ss << ".hot-"<<  couplingId; 
  fullFileName = ss.str();
}


void TopologyFile::initialize(const TreeAln& traln, nat someId, bool isDryRun)  
{    
  rejectIfExists(fullFileName);   

  if(isDryRun)
    return; 

  std::ofstream fh(fullFileName,std::fstream::out );  // std::fstream::app    
    
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

void TopologyFile::sample(const TreeAln &traln, nat gen,  AbstractParameter* param)  
{    
  auto&& fh = std::ofstream(fullFileName,std::fstream::app|std::fstream::out ); 
  auto tp = TreePrinter(true, false, false);
  
  auto treeString = tp.printTree(traln, param); 
  fh << "\ttree gen."<< gen
     << ".{" <<  param->getPartitions()  << "}"
     << " = [&U] " << treeString << std::endl; 

  fh.close(); 
}

void TopologyFile::finalize()  
{
  std::ofstream fh(fullFileName, std::fstream::app|std::fstream::out ); 
  fh << "end;" << std::endl; 
  fh.close(); 
}



void TopologyFile::regenerate(std::string workdir, std::string prevId, nat gen) 
{
  std::ofstream fh(fullFileName, std::fstream::out); 

  std::stringstream ss; 
  
  ss << OutputFile::getFileBaseName(workdir) << "_topologies." << prevId << "." << runid ; 
  if(hasManyTopoloFiles)
    ss << ".tree." << paramNum; 

  std::ifstream prevFile(ss.str()); 
  rejectIfNonExistant(ss.str()); 

  nat genFound = 0; 
  while(genFound < gen && not prevFile.eof())
    {
      std::string line; 
      getline(prevFile, line ); 
      
      // HACK 
      bool lineIsEmpty = line.compare("") == 0; 

      if(not lineIsEmpty)
	{
	  // find start
	  nat cnt = 0; 
	  while(line[cnt] == ' '  && cnt < line.size())
	    ++cnt;

	  std::stringstream ss; 
	  ss.str(line.substr(cnt+1)); 
	  // tout << "no white:>" << ss.str() << "<" << std::endl; 
	  std::string found;       
	  getline(ss, found, ' '); 

	  if(found.compare("tree") == 0)
	    {	  
	      getline(ss, found, ' ');       
	      genFound = std::stoi(found.substr(4)); 
	      // tout << "restoring file: at " << genFound << std::endl ; 
	    }

	  if(genFound < gen)
	    fh << line << std::endl ;       
	}
    }

  fh.close(); 
} 
