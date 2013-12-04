#include <iostream> 
#include <algorithm>
#include <cassert>
#include <cstring>
#include <vector>
#include <string>
#include <unistd.h>
#include <unordered_map> 
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iomanip> 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include "Arithmetics.hpp"
#include "common.h"

// #include "ParallelSetup.hpp"


void myExit(int code)
{
  exit(code); 
}

static bool isAAMod(const std::string& input )
{
  auto substr = input.substr(0,7) ; 
  return substr.compare("aaModel")  == 0 ; 
}


class Values
{
public: 
  size_t size() const {return std::max(values.size(), models.size()); }

  std::vector<double> values; 
  std::vector<std::string> models; 
}; 


std::unordered_map<std::string, Values> readFile(std::string file, nat burnin)
{
  auto result = std::unordered_map<std::string, Values>{}; 
  auto&& fh = std::ifstream(file); 

  auto line = std::string{}; 
  getline(fh, line);

  auto headers = std::vector<std::string>{} ; 
  getline(fh,line); 

  auto&& istr = std::istringstream(line); 
  auto elem = std::string{}; 
  while(getline(istr, elem, '\t'))
    headers.push_back(elem); 

  auto values = std::vector<Values>(); 

  for(nat i = 0; i < headers.size(); ++i)
    values.push_back(Values{}); 

  while(getline(fh, line))
    {
      nat ctr = 0; 
      auto &&istr = std::istringstream(line); 
      auto elem = std::string{}; 
      while(getline(istr, elem, '\t'))
	{
	  if(isAAMod(headers[ctr]))
	    {
	      values.at(ctr).models.push_back(elem); 
	    }
	  else 
	    {
	      auto &&elemHelper = std::istringstream(elem); 
	      double value = 0; 
	      elemHelper >> value; 
	      values.at(ctr).values.push_back(value);
	    }
	  ++ctr; 

	}
    }

  for(nat i = 0; i < headers.size() ; ++i)
    result[headers[i]] = values[i]; 

  // remove the burnin 
  if(burnin > 0 )
    {
      nat ctr = 0; 
      for(auto &elem : result)
	{
	  auto& values = std::get<1>(elem); 
	  
	  if(values.size() <= burnin)
	    {
	      std::cerr << "error: I merely have " << values.size() << " values and you are trying to skip " << burnin << " of them using the -b option" << std::endl; 
	      myExit(-1); 
	    }
	  
	  
	  std::cout << "querying " << headers[ctr] << std::endl; 
	  if(isAAMod(headers[ctr]))
	    values.models.erase(values.models.begin(), values.models.begin() + burnin);
	  else 
	    values.values.erase(values.values.begin(), values.values.begin() + burnin);

	  ++ctr ; 
	}
    } 

  return result; 
}


static void printUsage(std::ostream &out)
{
  out << "postProcParams is a utility that provides various statistics for sampled parameters.\n\n"; 
  out << "Usage: ./postProcParams -n id -f file[..] [-b burnin]\n\n"; 
  out << "    -n id             an id for the output file.\n"  ; 
  out << "    -b burnin         a number of samples to neglect (starting from the beginning of the file)\n"  ; 
  out << "    -f file[..]       one or many ExaBayes_parameters* files\n\n" ; 

}


std::tuple<std::string, std::vector<std::string>, nat> processCmdLine(int argc, char **argv)
{
  
  int c = 0; 
  auto id = std::string {}; 
  auto files = std::vector<std::string> {}; 
  auto burnin  = 0; 
  
  while( (c  = getopt(argc, argv, "hn:f:b:")) != EOF) 
    {
      switch(c)
	{
	case 'h': 
	  {
	    printUsage(std::cout); 
	    myExit(-1); 
	  }
	  break; 
	case 'n': 
	  {
	    id = std::string{optarg, optarg + strlen(optarg)}; 
	  }
	  break; 
	case 'f' : 
	  {
	    int index  = optind -1; 
	    while(index < argc)
	      {
		auto next = std::string{argv[index]}; 
		++index; 
		if(next[0] != '-') 
		  files.push_back(next);
		else 
		  {
		    optind =  index -1 ; 
		    break; 
		  }
	      }
	  }
	  break; 
	case 'b': 
	  {
	    auto &&ss = std::stringstream{optarg}; 
	    ss >> burnin; 
	  }
	  break; 
	default : 
	  {
	    std::cerr << "error: unknown argument "  << optarg << std::endl; 
	    myExit(-1); 
	  }
	}
    }

  if(files.size()  == 0)
    {
      std::cerr << "Error: please specify one or many files via -f." << std::endl; 
      myExit(-1); 
    }
  if(id.compare("") == 0)
    {
      std::cerr << "Error: please specify a run id via -n. " << std::endl; 
      myExit(-1); 
    }

  return std::make_tuple(id, files, burnin); 
}



int main(int argc, char **argv)
{
  if(argc < 2 )
    {
      printUsage(std::cout);
      myExit(-1); 
    }  

  auto files = std::vector<std::string>{}; 
  nat burnin = 0; 
  auto id = std::string{}; 
  std::tie(id, files, burnin) = processCmdLine(argc, argv);

  // open output file 
  auto && ss = std::stringstream{};
  ss <<  PROGRAM_NAME << "_parameterStatistics."  << id  ; 
  auto outputFileName = ss.str(); 
  if(std::ifstream(outputFileName))
    {
      std::cerr << "Error: output file >" << outputFileName << "< already exists. Please or move." << std::endl; 
      myExit(-1); 
    }
  
  auto && out = std::ofstream(outputFileName); 
  auto headerToValues = std::vector<std::unordered_map<std::string,Values>>{};
  for(auto &file : files )
    {
      std::cout << "reading file " << file  << std::endl; 
      headerToValues.push_back(readFile(file, burnin));
    }

  auto headers = std::vector<std::string>{}; 
  bool haveAA = false;
  for (auto elem :  headerToValues[0]) 
    {
      haveAA |= isAAMod(elem.first); 
      headers.push_back(elem.first); 
    }

  auto elemsToIgnore = std::vector<std::string>
    {
      "Gen", 
      // "LnPr",
      // "LnL"
    }; 
  

  nat numCols = 0; 
  out << "paramName\tmean\tsd\tperc5\tperc25\tmedian\tperc75\tper95\tESS" ;
  numCols += 9; 

  if(files.size() > 1)
    {
      out << "\tprsf"; 
      ++numCols; 
    }

  if(haveAA)
    {
      out << "\tdiscreteDistribution"; 
      ++numCols; 
    }

  out << std::endl; 
  
  
  for(auto header : headers)
    {
      bool isAA = isAAMod(header); 

      bool isGood = true; 
      for(auto elem : elemsToIgnore)
	isGood &= elem.compare(header) != 0 ; 
      if(not isGood)
	continue;  
      
      auto valuesConcat = Values{}; 
      auto relevant = std::vector<Values>{}; 
      for(auto headerToSomeValues : headerToValues)
	{
	  relevant.push_back(headerToSomeValues[header]); 
	  auto someVals = headerToSomeValues[header]; 

	  if(isAA)
	    {
	      valuesConcat.models.reserve(valuesConcat.size() + someVals.size()); 
	      valuesConcat.models.insert(end(valuesConcat.models), begin(someVals.models),end(someVals.models));
	    }
	  else 
	    {
	      valuesConcat.values.reserve(valuesConcat.size() + someVals.size()); 
	      valuesConcat.values.insert(end(valuesConcat.values), begin(someVals.values),end(someVals.values));
	    }
	}

      out << MAX_SCI_PRECISION ; 

      if(isAA)
	{
	  out << header; 
	  for(nat i = 0; i < numCols-1; ++i)
	    out << "\tNA" ;
	  
	  auto frequencyHash = std::unordered_map<std::string, double>{}; 
	  auto total = double{0}; 
	  for(auto &v : valuesConcat.models )
	    {
	      ++frequencyHash[v]; 
	      ++total; 
	    }
	  
	  for(auto &elem : frequencyHash)
	    elem.second /= total; 

	  bool isFirst = true; 
	  for(auto &elem : frequencyHash)
	    {
	      out << MAX_SCI_PRECISION << ( isFirst ? "" : ","  ) << elem.first << ":" << elem.second; 
	      isFirst = false; 
	    }

	  out << std::endl; 
	}
      else 			// print values 
	{
	  auto relVals = std::vector<std::vector<double>>{};
	  for(auto r : relevant)
	    relVals.push_back(r.values); 
      
	  auto prsf = Arithmetics::PRSF(relVals); 
	  auto sd = sqrt(Arithmetics::getVariance(valuesConcat.values));
	  auto perc95  = Arithmetics::getPercentile(.95, valuesConcat.values);
	  auto perc5 = Arithmetics::getPercentile(.5, valuesConcat.values); 
	  auto perc50 = Arithmetics::getPercentile(.50, valuesConcat.values); 
	  auto perc25 = Arithmetics::getPercentile(.25, valuesConcat.values); 
	  auto perc75 = Arithmetics::getPercentile(.75, valuesConcat.values); 
	  auto ess = Arithmetics::getEffectiveSamplingSize(valuesConcat.values); 
	  auto mean = Arithmetics::getMean(valuesConcat.values); 
	  out << header
	      << "\t" << mean
	      << "\t" << sd
	      << "\t" << perc5
	      << "\t" << perc25
	      << "\t" << perc50
	      << "\t" << perc75
	      << "\t" << perc95
	      << "\t" << ess ; 
      
	  if(relevant.size()  > 1)
	    out << "\t" << prsf; 

	  if(haveAA)
	    out << "\tNA" ; 

	  out << std::endl; 
	} 
    }

  std::cout << "successfully printed statistics to file >" <<  outputFileName << "<" << std::endl; 
  
  return 0; 
}
