#include "ParameterFile.hpp"
#include "GlobalVariables.hpp"

#include "ParallelSetup.hpp"

#include <cassert>


ParameterFile::ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId)
  : runid(runid)
  , couplingId(couplingId)
{
  std::stringstream ss ; 
  // TODO portability 
  ss << workdir << (workdir.compare("") == 0 ? "" : "/")  << PROGRAM_NAME << "_parameters." << runname << "." << runid ; 
  if(couplingId != 0 )
    ss << ".hot-"<<  couplingId; 
  fullFileName = ss.str();

  rejectIfExists(fullFileName); 
  tout << "initialized parameter file >" << fullFileName << "<" << std::endl; 

  std::ofstream fh(fullFileName); 
  fh << "" ; 
  fh.close(); 
}


void ParameterFile::initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId ) const 
{        
  std::ofstream fh(fullFileName,std::fstream::out);  // std::fstream::app|

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
  std::ofstream fh(fullFileName, std::fstream::app|std::fstream::out); 
    
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
  std::ofstream fh(fullFileName, std::fstream::out); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_parameters." << prevId << "." << runid ; 
  std::ifstream prevFile(ss.str()); 
  rejectIfNonExistant(ss.str()); 

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
