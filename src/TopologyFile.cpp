#include "TopologyFile.hpp"
#include "GlobalVariables.hpp"

#include <cassert>


TopologyFile::TopologyFile(std::string workdir, std::string runname, nat runid, nat couplingId) 
  : runid(runid)
  , couplingId(couplingId)
{ 
  std::stringstream ss ; 
  // TODO portability 
  ss << workdir <<  ( workdir.compare("") == 0  ? "" :  "/" )  << "ExaBayes_topologies." << runname << "." << runid ; 
  if(couplingId != 0 )
    ss << ".hot-"<<  couplingId; 
  fullFilename = ss.str();

  if( std::ifstream (fullFilename) ) 
    {
      std::cerr << std::endl <<  "File " << fullFilename << " already exists (probably \n"
		<< "from previous run). Please choose a new run-id or remove previous output files. " << std::endl; 
      exit(0); 
    }

  std::ofstream fh (fullFilename); 
  fh << "" ; 
  fh.close(); 
}


void TopologyFile::initialize(const TreeAln& traln, nat someId) const 
{    
  std::ofstream fh(fullFilename,std::fstream::out );  // std::fstream::app    
    
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




void TopologyFile::sample(const TreeAln &traln, nat gen) const 
{    
  std::ofstream fh(fullFilename,std::fstream::app|std::fstream::out ); 
  TreePrinter tp(true, false, false);
  std::string treeString = tp.printTree(traln);
  fh << "\ttree gen." << gen << " = [&U] " << treeString << std::endl; 
  fh.close(); 
}

void TopologyFile::finalize() const 
{
  std::ofstream fh(fullFilename, std::fstream::app|std::fstream::out ); 
  fh << "end;" << std::endl; 
  fh.close(); 
}



void TopologyFile::regenerate(std::string prevId, nat gen) 
{
  std::ofstream fh(fullFilename, std::fstream::out); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_topologies." << prevId << "." << runid ; 
  std::ifstream prevFile(ss.str()); 
  
  if(not prevFile)
    {
      std::cerr << "Error: could not find the topology file from previous run. \n"
		<< "The assumed name of this file was >" <<  ss.str() << "<. Aborting." << std::endl; 
      exit(0);
    }

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