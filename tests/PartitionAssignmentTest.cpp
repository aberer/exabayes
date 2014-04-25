#include "PartitionAssignment.hpp"
#include <algorithm>


class PartitionAssignmentTest : public testing::Test
{
public: 
  void SetUp()
  {
    auto states = std::vector<nat>
      {
	4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
      }; 
    
    auto length = std::vector<nat>
      {
	294, 469, 920,  69, 227,  63, 169, 258,  45, 278,  36, 466, 242, 153, 303, 241, 253, 333,  56, 157
      }; 

    nat numTax = 10; 
    for(nat i = 0; i< states.size();  ++i)
      {
	_partitions.emplace_back(numTax, "", PLL_DNA_DATA,states.at(i), 4, false );
	_partitions.back().setLower(0);
	_partitions.back().setUpper(length.at(i));
      }
  }; 

  void printStatistics( const PartitionAssignment &pa )
  {
    auto ass = pa.getAssignment();
    tout << "\nRESULT:\n"; 
    for(auto &a : ass)
      tout << std::get<1>(a) << std::endl;

    tout << "================================================================"<< std::endl ;
    auto numParts = pa.getNumPartPerProcess();
    auto numPats = pa.getSitesPerProcess();
    auto sum = std::accumulate(begin(numParts), end(numParts), 0u);
    for(nat i =0 ; i < pa.getNumProc(); ++i)
      tout << i << "\t" << numParts[i] << "\t" << numPats[i] <<  std::endl; 
    tout <<  " => "  << sum << " assignments for " << _partitions.size() << " partitions" << std::endl; 
  }


  void checkIfAllAssigned( const PartitionAssignment& pa)
  {
    nat totalSites = 0; 
    for(auto &p : _partitions)
      totalSites += p.getUpper() - p.getLower();

    nat check = 0; 
    auto ass = pa.getAssignment(); 
    for(auto &elem : ass)
      {
	auto assignment = elem.second;
	check += assignment.width; 
      }
    ASSERT_EQ(totalSites, check);
  }
  
  std::vector<Partition> _partitions; 
}; 


TEST_F(PartitionAssignmentTest, generic )
{
  for(nat i = 1; i < 2 ; ++i) // 100
    {
      auto pa =  PartitionAssignment(i); 
      pa.assign(_partitions);
      // printStatistics(pa);
      checkIfAllAssigned(pa); 
      // tout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl; 
    }
}




TEST(PartitionAssignment, tooFewSites)
{
  auto _partitions =  std::vector<Partition>() ; 
  _partitions.emplace_back(10, "", PLL_DNA_DATA,4, 4, false );
  _partitions.back().setLower(0);
  _partitions.back().setUpper(10); 
  
  auto pa = PartitionAssignment(12); 
  pa.assign(_partitions); 
}


