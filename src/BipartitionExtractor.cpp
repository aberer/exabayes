#include "BipartitionExtractor.hpp"
#include "ParallelSetup.hpp"

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
      ParallelSetup::genericExit(-1); 
    }
}


// does not include length of tips currently   
BipartitionExtractor::BipartitionExtractor(std::vector<std::string> files, bool extractToOneHash)
  : TreeProcessor(files)
  , _extractToOneHash(extractToOneHash)
{
  assert(taxa.size( ) != 0); 
}


BipartitionExtractor::BipartitionExtractor( BipartitionExtractor&& rhs) 
  : TreeProcessor(std::move(rhs))
  , _bipHashes(std::move(rhs._bipHashes))
  , _uniqueBips(std::move(rhs._uniqueBips))
  , _extractToOneHash(std::move(_extractToOneHash))
{
  // should work, but am skeptic
  assert(0); 
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
      ParallelSetup::genericExit(-1); 
    }

  result -= 2;
  return result; 
}


template<bool readBl>
void BipartitionExtractor::extractBips(nat burnin)
{
  int ctr = 0; 
  for (auto filename : fns)
    {
      auto &&ifh = std::ifstream(filename); 

      nat end = getNumTreesInFile(filename); 

      if( not _extractToOneHash || _bipHashes.size() == 0)
	_bipHashes.emplace_back(tralnPtr->getNumberOfTaxa());
      auto &bipHash = _bipHashes.back(); 

      // skip the burnin 
      nat i = 0; 
      for( ; i < burnin; ++i)
	skipTree(ifh); 

      for( ; i < end; ++i)
	{
	  nextTree<readBl>(ifh);
	  bipHash.addTree(*tralnPtr,true, true);
	}
      
      ++ctr; 
    }

  extractUniqueBipartitions();
}


void BipartitionExtractor::extractUniqueBipartitions()
{
  auto set = std::unordered_set<Bipartition>();

  for(auto &bipHash : _bipHashes)
    for(auto &bip : bipHash)
      set.insert(bip.first); 
  
  _uniqueBips.clear();
  nat ctr = 0; 
  for(auto &bip :set)
    {
      _uniqueBips[bip] = ctr; 
      ++ctr; 
    }
}


void BipartitionExtractor::printBipartitionStatistics(std::string id) const 
{
  assert(_uniqueBips.size() != 0); 

  auto && ss = std::stringstream{}; 
  ss << PROGRAM_NAME << "_bipartitionStatistics." << id ; 
  rejectIfExists(ss.str()); 

  auto && freqFile = std::ofstream (ss.str()); 
  freqFile << "id\t";

  freqFile
    << "freq"
    << "\tbl.mean"
    << "\tbl.sd"
    << "\tbl.ESS"
    << "\tbl.perc5"
    << "\tbl.perc25"
    << "\tbl.perc50"
    << "\tbl.perc75"
    << "\tbl.perc95"
    ; 

  if(fns.size() > 1)
    freqFile << "\tbl.prsf" ; 
  freqFile << std::endl; 

  for(auto &bipElem: _uniqueBips)
    {
      auto allBls = std::vector<std::vector<double>>{};      
      for(auto &bipHash : _bipHashes )
	allBls.push_back(bipHash.getBranchLengths(bipElem.first)); 

      auto allBlsConcat = std::vector<double>{}; 
      for(auto &allBl : allBls)
	allBlsConcat.insert(end(allBlsConcat), begin(allBl), end(allBl));
      
      assert(allBls.size() >  0); 

      auto num = allBlsConcat.size() ; 
      auto mean = Arithmetics::getMean(allBlsConcat) ; 
      auto sd = sqrt(Arithmetics::getVariance(allBlsConcat) ); 
      auto ess = Arithmetics::getEffectiveSamplingSize(allBlsConcat) ; 
      auto perc5= Arithmetics::getPercentile(0.05,allBlsConcat) ; 
      auto perc25 = Arithmetics::getPercentile(0.25,allBlsConcat)  ; 
      auto perc50 = Arithmetics::getPercentile(0.50,allBlsConcat) ; 
      auto perc75 = Arithmetics::getPercentile(0.75,allBlsConcat) ; 
      auto perc95 = Arithmetics::getPercentile(0.95,allBlsConcat) ; 

      freqFile << bipElem.second // the id 
	       << "\t" << num 
	       << "\t" << mean
	       << "\t" << sd 
	       << "\t" << ess 
	       << "\t" << perc5 
	       << "\t" << perc25
	       << "\t" << perc50 
	       << "\t" << perc75
	       << "\t" << perc95
	; 
	
      if(fns.size() > 1)
	freqFile << "\t" << Arithmetics::PRSF(allBls); 

      freqFile << std::endl; 
    }

  
  std::cout << "printed bipartition statistics to file " << ss.str() << std::endl; 
}


void BipartitionExtractor::printBranchLengths(std::string id) const 
{
  assert(_uniqueBips.size() != 0); 
  
  auto &&ss = std::stringstream {}; 
  ss << PROGRAM_NAME << "_bipartitionBranchLengths."  << id; 
  rejectIfExists(ss.str()); 

  std::ofstream blFile(ss.str());

  blFile << MAX_SCI_PRECISION; 
  blFile  << "bipId\tfileId\tlength"  << std::endl;
  nat ctr = 0; 
  for(auto &bipHash : _bipHashes)
    {
      for(auto &bip : bipHash)
	for(auto length : bipHash.getBranchLengths(bip.first) )
	  blFile << _uniqueBips.at(bip.first) << "\t"  << ctr  << "\t" << length << std::endl; 
      ++ctr ; 
    }

  std::cout << "printed branch lengths associated with bipartitions to " << ss.str() << std::endl; 
}

void BipartitionExtractor::printFileNames(std::string id) const 
{
  assert(_uniqueBips.size() != 0); 

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
  assert(_uniqueBips.size() != 0); 

  std::stringstream ss; 
  ss << PROGRAM_NAME << "_bipartitions." << id; 
  rejectIfExists(ss.str()); 

  std::ofstream out(ss.str());

  for(auto &bipElem :  _uniqueBips)
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



void BipartitionExtractor::buildTreeRecursive(nat currentId, const std::vector<std::vector<nat> >& directSubBips, 
					      const std::vector<Bipartition> &bips, std::stringstream &result, bool printSupport) const
{
  nat totalTrees = 0;
  for(auto &bipHash : _bipHashes)
    totalTrees += bipHash.getTreesAdded(); 

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
      double mySupport = 0; 
      for(auto &bipHash : _bipHashes)
	mySupport += bipHash.getPresence(curBip).count() ; 
      mySupport /= double(totalTrees); 

      assert(mySupport > 0. && mySupport <= 1.); 

      result << ")" ; 
      if(printSupport)
      	result << ":[" << std::fixed << std::setprecision(3) << mySupport << "]"; 
    }
  else 
    result << ");"; 
}


template void BipartitionExtractor::extractBips<true>(nat burnin); 
template void BipartitionExtractor::extractBips<false>(nat burnin); 

