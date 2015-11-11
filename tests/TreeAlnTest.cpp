#include "model/TreeAln.hpp" 
#include "system/ByteFile.hpp"


TEST(tree, something)
{
  nat numTax = 10; 
  nat numPart = 3; 

  auto&& traln = TreeAln(numTax, false); 
  traln.createCaterpillar();
  
  
  auto bl = BranchLengthResource{}; 
  bl.initialize(numTax, numPart); 
  
  traln.setBranchLengthResource(bl);

  auto&& traln2 = TreeAln(numTax, false); 
  traln2 = traln; 
  
}


TEST(tree, assignment)
{
  // auto&& bf = ByteFile("/home/aberer/proj/exa-bayes/data/tiny/aln.binary"); 
  

}
