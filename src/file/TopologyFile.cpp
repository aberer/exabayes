#include "TopologyFile.hpp"
#include "system/GlobalVariables.hpp"
#include "system/extensions.hpp"
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

  auto taxa =  traln.getTaxa(); 
  nat ctr = 1; 
  for(auto taxon : taxa)
    {
      auto isLast = taxon.compare( taxa.back() )  == 0 ; 
      fh << "\t" <<  ctr << "\t"
  	 << taxon  << (isLast ? ";" : ",")<< std::endl; 
      ++ctr; 
    }

  fh.close(); 
}


std::streamoff TopologyFile::getPosBeforeEnd() const 
{
  auto &&in = std::ifstream{fullFileName};
  // in.sync();
  
  int c = 0; 
  int pos = 0; 
  nat newLinesEncountered = 0; 
  
  while(newLinesEncountered < 2 && c != EOF)
    {
      --pos; 
      in.seekg(pos, std::ios::end); 
      c = in.peek();
      if(c == '\n')
	++newLinesEncountered; 
    }
  
  ++pos; 
  in.seekg(pos, std::ios::end); 
  auto str = std::string(); 
  getline(in, str); 

  if(str.compare("end;") == 0)
    return std::streamoff{pos}; 
  else 
    return 0; 
}


void TopologyFile::sample(const TreeAln &traln, nat gen,  AbstractParameter* param)  
{    
  auto pos = getPosBeforeEnd(); 

  auto&& fh = std::fstream(fullFileName, std::ios::in | std::ios::out);
  fh.seekp(pos, std::ios::end); 

  auto tp = TreePrinter( param->getCategory() == Category::BRANCH_LENGTHS , false, false);

  auto treeString = tp.printTree(traln, param);
  fh << "\ttree gen."<< gen
     << ".{"; 

  formatRange(fh, param->getPartitions()); 

  fh << "}"
     << " = [&U] " << treeString << std::endl; 

  fh << "end;" << std::endl; 

}


void TopologyFile::regenerate(std::string workdir, std::string prevId, nat gen) 
{
  auto&& fh = std::ofstream(fullFileName, std::fstream::out); 

  auto&& ss = std::stringstream {}; 
  
  ss << OutputFile::getFileBaseName(workdir) << "_topologies." << prevId << "." << runid ; 
  if(hasManyTopoloFiles)
    ss << ".tree." << paramNum; 

  auto&& prevFile = std::ifstream(ss.str()); 
  rejectIfNonExistant(ss.str()); 

  nat genFound = 0; 
  while(genFound < gen &&  prevFile)
    {
      auto line = std::string{}; 
      getline(prevFile, line ); 
      
      // HACK 
      bool lineIsEmpty = line.compare("") == 0; 

      if(not lineIsEmpty)
	{
	  // find start
	  nat cnt = 0; 
	  while(line[cnt] == ' '  && cnt < line.size())
	    ++cnt;

	  auto&& ss = std::stringstream{} ; 
	  ss.str(line.substr(cnt+1)); 
	  // tout << "no white:>" << ss.str() << "<" << std::endl; 
	  auto found = std::string{} ;       
	  getline(ss, found, ' '); 

	  if(found.compare("tree") == 0)
	    {	  
	      getline(ss, found, ' ');       
	      genFound = std::stoi(found.substr(4)); 
	      // tout << "restoring file: at " << genFound << std::endl ; 
	    }

	  if(genFound <= gen)
	    fh << line << std::endl ;       
	}
    }
} 
