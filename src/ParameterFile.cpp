#include "ParameterFile.hpp"
#include "GlobalVariables.hpp"

#include <cassert>


ParameterFile::ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId)
  : runid(runid)
  , couplingId(couplingId)
{
  std::stringstream ss ; 
  // TODO portability 
  ss << workdir << (workdir.compare("") == 0 ? "" : "/")  << "ExaBayes_parameters." << runname << "." << runid ; 
  if(couplingId != 0 )
    ss << ".hot-"<<  couplingId; 
  fullFilename = ss.str();

  if( std::ifstream(fullFilename) ) 
    {
      std::cerr << std::endl <<  "File " << fullFilename << " already exists (probably \n"
		<< "from previous run). Please choose a new run-id or remove previous output files. " << std::endl; 
      exit(0); 
    }

  std::ofstream fh(fullFilename); 
  fh << "" ; 
  fh.close(); 
}


void ParameterFile::initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId ) const 
{        
  std::ofstream fh(fullFilename,std::fstream::out);  // std::fstream::app|

  fh << "[ID: " << someId << "]" << std::endl; 

  fh << "Gen\t";
  fh << "LnPr\t"; 
  fh << "LnL\t" ; 
  fh << "TL\t" ; 

  bool isFirst = true; 
  for(auto &p : parameters)
    {
      if(p->isPrintToParamFile())
	{
	  if(isFirst) 
	    isFirst = false; 
	  else 
	    fh << "\t" ; 
	  p->printAllComponentNames(fh, traln); 
	}
    }

  fh << std::endl; 

  fh.close(); 
}


void ParameterFile::sample(const TreeAln &traln, const std::vector<AbstractParameter*> parameters, nat gen, double lnPr) const 
{
  std::ofstream fh(fullFilename, std::fstream::app|std::fstream::out); 
    
  fh << gen << "\t"; 
  fh << std::setprecision(std::numeric_limits<double>::digits10)  << std::scientific; 
  fh << lnPr << "\t"; 
  fh << traln.getTr()->likelihood << "\t" ; 
  fh << Branch(0,0,traln.getTreeLengthExpensive()).getInterpretedLength(traln) << "\t"; 

  bool isFirst = true; 
  for(auto &p : parameters)
    {
      if(p->isPrintToParamFile())
	{
	  if(isFirst)
	    isFirst = false; 
	  else 
	    fh << "\t"; 
	  p->printSample(fh,traln);
	}
    }
  fh << std::endl; 
    
  fh.close(); 
}


void ParameterFile::regenerate(std::string prevId, nat gen) 
{
  std::ofstream fh(fullFilename, std::fstream::out); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_parameters." << prevId << "." << runid ; 
  std::ifstream prevFile(ss.str()); 

  if(not prevFile)
    {
      std::cerr << "Error: could not find the parameter file from previous run. \n"
		<< "The assumed name of this file was >" <<  ss.str() << "<. Aborting." << std::endl; 
      exit(0);
    }

  nat genFound = 0; 
  nat lineCtr = 0; 
  while(genFound < gen && not prevFile.eof())
    {
      std::string line; 
      getline(prevFile, line); 
      
      // hack =/ 
      bool lineIsEmpty = line.compare("") == 0 ; 
      
      if(lineCtr > 1 && not lineIsEmpty)
	{
	  std::stringstream ss; 
	  ss.str(line); 
	  // tout << "line=" << line << std::endl; 
	  std::string part; 
	  getline( ss ,part, '\t'  );       
	  genFound = std::stoi(part); 
	}

      if(genFound < gen && not lineIsEmpty)
	fh << line << std::endl; 
      ++lineCtr; 
    }

  fh.close();
} 
