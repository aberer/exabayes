#include "TreeAln.hpp" 


TEST(tree, something)
{
  nat numTax = 10; 
  nat numPart = 3; 

  auto traln = TreeAln(numTax); 
  traln.createCaterpillar();
  
  
  auto bl = BranchLengthResource{}; 
  bl.initialize(numTax, numPart); 
  
  traln.setBranchLengthResource(bl);

  auto traln2 = TreeAln(numTax); 
  traln2 = traln; 
  
}
