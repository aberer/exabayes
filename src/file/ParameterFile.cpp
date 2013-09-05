#include "ParameterFile.hpp"
#include "GlobalVariables.hpp"
#include "Category.hpp"

#include "ParallelSetup.hpp"

#include <cassert>


ParameterFile::ParameterFile(std::string workdir, std::string runname, nat runid, nat couplingId)
  : runid(runid)
  , couplingId(couplingId)
{
  std::stringstream ss ; 

  ss << OutputFile::getFileBaseName(workdir) << "_parameters." << runname << "." << runid ; 
  if(couplingId != 0 )
    ss << ".hot-"<<  couplingId; 
  fullFileName = ss.str();

}


void ParameterFile::initialize(const TreeAln& traln, std::vector<AbstractParameter*> parameters,  nat someId )  
{        
  rejectIfExists(fullFileName); 
  tout << "initialized parameter file >" << fullFileName << "<" << std::endl; 

  std::ofstream fh(fullFileName,std::fstream::out);  // std::fstream::app|

  fh << "[ID: " << someId << "]" << std::endl; 

  fh << "Gen\t";
  fh << "LnPr\t"; 
  fh << "LnL\t" ; 

  std::vector<AbstractParameter*>  blParams; 
  for(auto &param : parameters)
    if(param->getCategory() == Category::BRANCH_LENGTHS)
      blParams.push_back(param); 

  for(auto &param : blParams)    
    fh << "TL{" << param->getPartitions()   << "}\t"; 

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


void ParameterFile::sample(const TreeAln &traln, const std::vector<AbstractParameter*> parameters, nat gen, double lnPr)  
{
  std::ofstream fh(fullFileName, std::fstream::app|std::fstream::out); 
    
  fh << gen << "\t"; 
  fh << MAX_SCI_PRECISION; 
  fh << lnPr << "\t"; 
  fh << traln.getTr()->likelihood << "\t" ; 

  std::vector<AbstractParameter*> blParams ; 
  for(auto &p : parameters)
    {
      if(p->getCategory() == Category::BRANCH_LENGTHS)
	blParams.push_back(p);
    }
  // should be possible 
  // copy_if(  parameters.begin(), parameters.end(), blParams.begin(), [=](AbstractParameter* p){ return p->getCategory() == Category::BRANCH_LENGTHS;  }); 

  // fh << Branch(0,0,traln.getTreeLengthExpensive()).getInterpretedLength(traln) << "\t"; 

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


void ParameterFile::regenerate(std::string workdir, std::string prevId, nat gen) 
{
  std::ofstream fh(fullFileName, std::fstream::out); 

  std::stringstream ss; 
  ss << OutputFile::getFileBaseName(workdir) << "_parameters." << prevId << "." << runid ; 
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
