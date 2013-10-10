#include "BipartitionExtractor.hpp"

#include <iostream>
#include <unordered_set>
#include <cassert>
#include <fstream>

#include "Arithmetics.hpp" 
#include "BipartitionHash.hpp"

#include "TreePrinter.hpp"
#include "parameters/BranchLengthsParameter.hpp"


static void rejectIfExists(std::string filename)
{
  if( std::ifstream(filename) ) 
    {
      std::cerr << std::endl <<  "File " << filename << " already exists (probably \n"
		<< "from previous run). Please choose a new run-id or remove previous output files. " << std::endl; 
      exit(-1); 
    }
}


// does not include length of tips currently   
BipartitionExtractor::BipartitionExtractor(std::vector<std::string> files)
  : TreeProcessor(files)
{
  assert(taxa.size( ) != 0); 
}


BipartitionExtractor::BipartitionExtractor( BipartitionExtractor&& rhs) 
  : TreeProcessor(std::move(rhs))
  , bipHashes(std::move(rhs.bipHashes))
  , uniqueBips(std::move(rhs.uniqueBips))
{
}


BipartitionExtractor& BipartitionExtractor::operator=(BipartitionExtractor rhs)
{
  if(this == &rhs )
    return *this; 
  else 
    {
      assert(0); 
    }
}


nat BipartitionExtractor::getNumTreesInFile(std::string file) const 
{
  // VERY UNSAFE  

  nat result = 0; 
  
  std::ifstream fh(file); 
  auto line = std::string();
  while(getline(fh, line))
    {
      if(line.compare("end;") == 0 )
	break;

      if(line.find(";") != std::string::npos) 
	++result; 
    }

  if(result < 1)
    {
      std::cout << "Header of file " << file << " possibly broken. Expected list of taxa terminated by semi-colon." << std::endl; 
      exit(-1); 
    }

  result -= 2;
  return result; 
}


void BipartitionExtractor::extractBipsNew()
{
  int ctr = 0; 
  for (auto filename : fns)
    {
      auto &&ifh = std::ifstream(filename); 

      nat end = getNumTreesInFile(filename); 

      auto bipHash = BipartitionHashNew(traln->getNumberOfTaxa()); 

      for(nat i = 0 ; i < end; ++i)
	{
	  nextTree(ifh);
	  bipHash.addTree(*traln,true);
	}
      bipHashes.push_back(bipHash);
 
      ++ctr; 
    }

  extractUniqueBipartitions();
}


void BipartitionExtractor::extractUniqueBipartitions()
{
  auto set = std::unordered_set<Bipartition>();

  for(auto &bipHash : bipHashes)
    for(auto &bip : bipHash)
      set.insert(bip.first); 
  
  uniqueBips.clear();
  nat ctr = 0; 
  for(auto &bip :set)
    {
      uniqueBips[bip] = ctr; 
      ++ctr; 
    }
}



void BipartitionExtractor::printBipartitionStatistics(std::string id) const 
{
  assert(uniqueBips.size() != 0); 

  auto && ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_bipartitionStatistics." << id ; 
  rejectIfExists(ss.str()); 

  auto && freqFile = std::ofstream (ss.str()); 
  freqFile << "id\t";
  for(auto &fn : fns)
    freqFile
      << "freq." << fn 
      << "\tbl.mean." << fn
      << "\tbl.sd." << fn
      << "\tbl.ESS."<< fn
      << "\t" ; 

  if(fns.size() > 1)
    freqFile << "prsf" ; 
  freqFile << std::endl; 

  auto allBlSets = std::vector<std::vector<double>>{};

  for(auto &bipElem: uniqueBips)
    {
      freqFile << bipElem.second << "\t" ; 
      for(auto &bipHash : bipHashes )
	{
	  auto bls = bipHash.getBranchLengths(bipElem.first); 
	  allBlSets.push_back(bls); 

	  auto mean = Arithmetics::getMean(bls);
	  auto var = Arithmetics::getVariance(bls); 
	  auto ess = Arithmetics::getEffectiveSamplingSize(bls);

	  freqFile << bipHash.getPresence(bipElem.first).count()
		   << "\t"  << mean
		   << "\t" << sqrt(var)
		   << "\t" << ess 
		   << "\t"
	    ; 
	}
      
      if(fns.size() > 1)
	freqFile << Arithmetics::PRSF(allBlSets); 

      freqFile << std::endl; 
    }

  
  std::cout << "printed bipartition statistics to file " << ss.str() << std::endl; 
}





void BipartitionExtractor::printBranchLengths(std::string id) const 
{
  assert(uniqueBips.size() != 0); 
  
  auto &&ss = std::stringstream {}; 
  ss << PROGRAM_NAME << "_bipartitionBranchLengths."  << id; 
  rejectIfExists(ss.str()); 

  std::ofstream blFile(ss.str());

  blFile << MAX_SCI_PRECISION; 
  blFile  << "bipId\tfileId\tlength"  << std::endl;
  nat ctr = 0; 
  for(auto &bipHash : bipHashes)
    {
      for(auto &bip : bipHash)
	for(auto length : bipHash.getBranchLengths(bip.first) )
	  blFile << uniqueBips.at(bip.first) << "\t"  << ctr  << "\t" << length << std::endl; 
      ++ctr ; 
    }

  std::cout << "printed branch lengths associated with bipartitions to " << ss.str() << std::endl; 
}

void BipartitionExtractor::printFileNames(std::string id) const 
{
  assert(uniqueBips.size() != 0); 

  auto &&ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_fileNames." << id; 
  rejectIfExists(ss.str()); 

  auto  &&ff = std::ofstream(ss.str());
  ff << "id\tfileName" << std::endl; 
  nat ctr = 0; 
  for(auto &fn : fns)
    {
      ff << ctr << "\t" << fn << std::endl; 
      ++ctr; 
    }

  std::cout << "printed file name identifiers (for future reference) to "<< ss.str() << std::endl; 
}


void BipartitionExtractor::printBipartitions(std::string id) const 
{
  assert(uniqueBips.size() != 0); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_bipartitions." << id; 
  rejectIfExists(ss.str()); 

  std::ofstream out(ss.str());

  for(auto &bipElem :  uniqueBips)
    {
      out << bipElem.second << "\t" ; 
      bipElem.first.printVerbose(out, taxa);
      out << std::endl;       
    }

  std::cout << "printed bipartitions and identifiers to file " << ss.str() << std::endl; 
}


std::string BipartitionExtractor::bipartitionsToTreeString(std::vector<Bipartition> bips, bool printSupport) const 
{
  // contains ids of direct children 
  auto directSubBips  = std::vector<std::vector<nat> >(bips.size()); 

  std::sort(bips.begin(), bips.end(),
	    [](const Bipartition& bipA, const Bipartition& bipB)
	    {
	      return bipA.count() < bipB.count() ; 
	    }
	    ); 

  auto toplevelBip = Bipartition();
  const auto& taxa = getTaxa();	
  toplevelBip.reserve(taxa.size()); 

  for(nat i = 0; i < taxa.size(); ++i)
    toplevelBip.set(i); 
  auto belowTop = std::vector<nat>{}; 

  // search the parent bipartition for each bipartition 
  for(nat i = 0; i < bips.size(); ++i)
    {
      bool parentFound = false; 
      const auto &bipA = bips[i]; // a <=> i
      for(nat j = i+1 ; j < bips.size() ;++j)
	{
	  const auto &bipB = bips[j];  // b <=> j 
	  if(bipA.isSubset(bipB))
	    {

	      directSubBips[j].push_back(i);
	      parentFound = true ; 
	      break; 
	    }
	}
      
      if( not parentFound )
	belowTop.push_back(i);
    }

  bips.push_back(toplevelBip); 
  directSubBips.push_back(belowTop); 

  std::stringstream ss; 
  buildTreeRecursive(bips.size()-1, directSubBips, bips, ss, printSupport); 
 
  return ss.str(); 
} 



void BipartitionExtractor::buildTreeRecursive(nat currentId, const std::vector<std::vector<nat> >& directSubBips, const std::vector<Bipartition> &bips, std::stringstream &result, bool printSupport) const
{
  assert(getBipartitionHashes().size() == 1); 
  nat totalTrees = getBipartitionHashes()[0].getTreesAdded(); 

  const auto& curBip = bips.at(currentId); 
  const auto& taxa = getTaxa();	// meh: efficiency 
  result << "("; 

  auto children = directSubBips.at(currentId); 
  auto taxaSeen = Bipartition{}; 
  taxaSeen.reserve(taxa.size());
  
  // print clades underneath this one  
  bool isFirst = true; 
  for(auto childId : children )
    {
      auto &childBip = bips.at(childId); 
      if( isFirst)
	isFirst = false; 
      else 
	result << ","; 
	
      buildTreeRecursive(childId, directSubBips, bips, result, printSupport); 
      taxaSeen = taxaSeen | childBip; 
    }
  
  if(children.size( ) != 0 && taxaSeen != curBip)
    result << ","; 

  // print remaining taxa in the bipartition that have not been
  // covered yet by sub-bipartitions (=> multifurcations)
  isFirst = true; 
  nat ctr = 0; 

  for(nat i = 0; i < curBip.getElemsReserved() && i < taxa.size() ;++i)
    {
      if( curBip.isSet(i) && not taxaSeen.isSet(i))
	{
	  if( isFirst)
	    isFirst = false; 
	  else 
	    result << ","; 
	    
	  result << taxa[i] ;
	  ++ctr; 
	}
    }

  if(currentId != bips.size()  -1 )
    {
      double mySupport = double(getBipartitionHashes()[0].getPresence(curBip).count()) / double(totalTrees); 
      result << ")" ; 
      if(printSupport)
      	result << ":[" << std::fixed << std::setprecision(3) << mySupport << "]"; 
    }
  else 
    result << ");"; 
}
