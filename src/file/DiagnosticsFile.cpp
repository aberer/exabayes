#include "DiagnosticsFile.hpp"

#include <sstream>
#include <unordered_map>


std::string DiagnosticsFile::createName(std::string runname, std::string workdir)
{
  std::stringstream ss; 
  ss << workdir << (workdir.compare("") == 0 ? "" : "/") 
     << PROGRAM_NAME << "_diagnostics." << runname ; 
  
  return ss.str();
}


void DiagnosticsFile::initialize(std::string workdir, std::string name, const std::vector<CoupledChains> &runs) 
{
  assert(not initialized); 
  initialized = true; 

  fullFileName = createName(name, workdir);     
  rejectIfExists(fullFileName); 

  // tout << "initialized diagnostics file >" << fullFileName << "<" << std::endl; 

  std::ofstream fh(fullFileName); 

  fh << "GEN"
     << "\tasdsf"; 
  
  for(auto& run : runs)
    {
      auto numElem = run.getChains().size(); 
      auto info = run.getSwapInfo(); 
      for(nat i = 0; i < numElem; ++i)
	{
	  for(nat j = i+1; j < numElem; ++j )
	    {
	      // auto sctr = info.getCounter(i,j); 
	      fh << "\tsw("  << i << "," << j << ")$run" << run.getRunid(); 
	    }
	}
    }
  
  for(auto& run : runs)
    {
      auto ps = run.getChains()[0].getProposalView(); 
      std::stringstream ss; 
      for(auto &p : ps)
	{
	  ss.str(""); 
	  p->printShort(ss); 
	  ss << "$run" << run.getRunid();  
	  names.push_back(ss.str()); 
	  fh << "\t" << ss.str(); 
	}        
    }

  fh << std::endl; 
  fh.close(); 
}


void DiagnosticsFile::printDiagnostics(nat gen, double asdsf, const std::vector<CoupledChains> &runs ) 
{
  std::ofstream fh(fullFileName, std::fstream::app); 

  fh << gen << "\t" << asdsf;   

  // print swapping info (if applicable)
  for(auto &run : runs)
    {
      auto m = run.getSwapInfo().getMatrix();
      for(auto &elem  : m)
	{
	  fh << "\t" << elem.getRatioInLastInterval()
	     << "," << elem.getRatioOverall() 
	     << "," << elem.getTotalSeen();
	}
    }

  // print acceptance rates  
  std::unordered_map<std::string,AbstractProposal*> name2proposal; 
  for(auto &run : runs)
    {
      for(auto &chain: run.getChains())
	{
	  if(chain.getCouplingId() != 0 ) // no hot chains
	    continue; 

	  auto &ps = chain.getProposalView(); 
	  for(auto &p : ps)
	    {
	      std::stringstream ss; 	      
	      p->printShort(ss);
	      ss << "$run"  << run.getRunid(); 

	      assert(name2proposal.find(ss.str()) == name2proposal.end()); 
	      name2proposal[ss.str()] = p; 
	    }
	}
    }

  for(auto &name : names) 
    {
      if(name2proposal.find(name) == name2proposal.end())
	{
	  tout << "could not find proposal " << name << std::endl; 
	  assert(0); 
	}
      auto& p = name2proposal[name]; 
      auto &sctr = p->getSCtr();
      fh << "\t" << 
	// sctr.getRatioInLast100()
	sctr.getRatioInLastInterval()
	 << "," << sctr.getRatioOverall() << "," << sctr.getTotalSeen() ; 
    }	  

  fh << std::endl; 
  fh.close(); 
}


void DiagnosticsFile::regenerate(std::string workdir, std::string nowId, std::string prevId, nat gen, nat numSwapEntries)
{    
  assert(not initialized); 
  initialized = true; 

  fullFileName = createName (nowId, workdir); 
  rejectIfExists(fullFileName); 
  std::ofstream fh(fullFileName);   
  std::string prevFileName = createName(prevId, workdir); 
  rejectIfNonExistant(prevFileName); 
  std::ifstream ifh(prevFileName); 

  nat genFound = 0; 
  nat lineCtr = 0; 
  bool firstLine = true;  

  while(genFound < gen && not ifh.eof())
    {
      std::string line ; 
      getline(ifh, line); 

      // special treatment
      if(firstLine)
	{
	  firstLine = false; 
	  std::stringstream ss; 
	  ss.str(line); 
	  
	  nat ctr = 0; 
	  std::string part; 
	  bool isEmpty = false; 
	  do 
	    {	      
	      auto& ret = getline(ss,part, '\t'); 
	      isEmpty = not ret; 
	      
	      // skipping over the swap matrix and additional info 
	      if( not (ctr < numSwapEntries + 2  ))
		{
		  // tout << "parsed " << part << std::endl; 
		  names.push_back(part); 
		}

	      ++ctr ; 
	    }while(part.compare("") !=  0 && not isEmpty); 	  
	}

      // TODO: all this comparing against empty lines only is needed,
      // because we could not have any data lines (and only the
      // header). Come up with a better solution. 
      if(lineCtr > 0 &&  not line.compare("") == 0
	 )
      	{
	  std::stringstream ss; 
	  ss.str(line); 
	  std::string part; 
	  getline(ss, part, '\t'); 
	  genFound = std::stoi(part); 
	}

      if(genFound < gen && not line.compare("") == 0)
	fh << line << std::endl; 
      ++lineCtr;
    }
  
  fh.close(); 
} 



