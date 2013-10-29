#include "ProtModelParameter.hpp" 
#include "ProtModel.hpp"


void ProtModelParameter::applyParameter(TreeAln& traln,  const ParameterContent &content) const 
{
  assert(content.protModel.size() == 1); 
  // auto raxmlValue = int(content.protModel); 
  // assert(0);
  // TODO  

  std::cerr << "\n\nTODO apply the aa parameter! (ProtModelParameter.cpp)\n\n" << std::endl; 
}
 
ParameterContent ProtModelParameter::extractParameter(const TreeAln &traln)  const  
{
  auto partition = traln.getPartition(partitions.at(0)); 
  auto model = ProtModel(partition->protModels); 
  auto result =   ParameterContent{} ; 
  result.protModel = {model}; 
  return result; 
}
   
void ProtModelParameter::printSample(std::ostream& fileHandle, const TreeAln &traln ) const 
{
  auto partition = traln.getPartition(partitions.at(0)); 
  auto name = ProtModelFun::getName(ProtModel(partition->protModels )); 
  fileHandle << name << "\t" ; 
}
 
void ProtModelParameter::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  
{
  fileHandle << "aaModel{" ; 
  bool isFirst = true; 
  for(auto &p : partitions)
    {
      fileHandle << (isFirst ? "" : ",") << p ; 
      isFirst = false; 
    }
  fileHandle << "}"; 
}

 
void ProtModelParameter::verifyContent(const TreeAln &traln, const ParameterContent &content) const  
{
  if(content.protModel.size() != 1)
    {
      tout << "incorrect number of models in parameter content. This is a programming error." << std::endl; 
      assert(0); 
    }
} 
