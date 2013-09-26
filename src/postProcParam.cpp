#include <iostream> 
#include <unordered_map> 
#include <fstream>
#include <algorithm>
#include <sstream>
#include <iomanip> 

#include "Arithmetics.hpp"
#include "common.h"


std::unordered_map<std::string, std::vector<double>> readFile(std::string file)
{
  auto result = std::unordered_map<std::string, std::vector<double>>{}; 
  std::ifstream fh(file); 

  auto line = std::string{}; 
  getline(fh, line);

  auto headers = std::vector<std::string>{} ; 
  getline(fh,line); 

  std::istringstream istr(line); 
  auto elem = std::string{}; 
  while(getline(istr, elem, '\t'))
    headers.push_back(elem); 

  auto values = std::vector<std::vector<double>>(); 

  for(nat i = 0; i < headers.size(); ++i)
    values.push_back(std::vector<double>{}); 

  while(getline(fh, line))
    {
      nat ctr = 0; 
      std::istringstream istr(line); 
      auto elem = std::string{}; 
      while(getline(istr, elem, '\t'))
	{
	  auto &&elemHelper = std::istringstream(elem); 
	  double value = 0; 
	  elemHelper >> value; 
	  values.at(ctr).push_back(value);
	  ++ctr; 
	}
    }

  // TODO watch out for amino acid stuff 

  for(nat i = 0; i < headers.size() ; ++i)
    result[headers[i]] = values[i]; 
    
  return result; 
}


int main(int argc, char **argv)
{
  if(argc < 2 )
    {
      std::cerr << "Usage: ./postProcParams file[..]\n"; 
      std::cerr << "\twhere file(s) can be multiple ExaBayes parameter files." << std::endl; 
      exit(-1); 
    }  

  auto files = std::vector<std::string>{}; 
  for(nat i = 1; i < nat(argc) ; ++i)
    files.push_back(std::string(argv[i]));

  auto headerToValues = std::vector<std::unordered_map<std::string,std::vector<double>>>{};
  for(auto &file : files )
    {
      std::cout << "reading file " << file  << std::endl; 
      headerToValues.push_back(readFile(file));
    }

  auto headers = std::vector<std::string>{}; 
  for (auto elem :  headerToValues[0]) 
    headers.push_back(elem.first); 

  auto elemsToIgnore = std::vector<std::string>
    {
      "Gen", 
      "LnPr",
      "LnL"
    }; 

  std::cout << "paramName\tmean\tsd\tperc5\tperc25\tmedian\tperc75\tper95\tprsf" << std::endl; 
  
  for(auto header : headers)
    {
      bool isGood = true; 
      for(auto elem : elemsToIgnore)
	isGood &= elem.compare(header) != 0 ; 
      if(not isGood)
	 continue;  
      
      auto valuesConcat = std::vector<double>{}; 
      auto relevant = std::vector<std::vector<double>>{}; 
      for(auto headerToSomeValues : headerToValues)
	{
	  relevant.push_back(headerToSomeValues[header]); 
	  auto someVals = headerToSomeValues[header]; 
	  valuesConcat.reserve(valuesConcat.size() + someVals.size()); 
	  valuesConcat.insert(valuesConcat.end(), someVals.begin(),someVals.end());
	}

      // std::cout << "concat: " ; 
      // for(auto &v : valuesConcat)
      // 	std::cout << v << ","; 
      // std::cout << std::endl; 

      std::cout << MAX_SCI_PRECISION ; 

      auto prsf = Arithmetics::PRSF(relevant); 
      auto sd = sqrt(Arithmetics::getVariance(valuesConcat));
      auto perc95  = Arithmetics::getPercentile(.95, valuesConcat);
      auto perc5 = Arithmetics::getPercentile(.5, valuesConcat); 
      auto perc50 = Arithmetics::getPercentile(.50, valuesConcat); 
      auto perc25 = Arithmetics::getPercentile(.25, valuesConcat); 
      auto perc75 = Arithmetics::getPercentile(.75, valuesConcat); 
      auto mean = Arithmetics::getMean(valuesConcat); 
      std::cout << header
		<< "\t" << mean
		<< "\t" << sd
		<< "\t" << perc5
		<< "\t" << perc25
		<< "\t" << perc50
		<< "\t" << perc75
		<< "\t" << perc95
		<< "\t" << prsf
		<< std::endl; 
    }
  
  return 0; 
}

