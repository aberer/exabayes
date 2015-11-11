#ifndef REAL_PARSER

/** 
    This file just forces the legacy parser into a c++ class 
*/ 

#include <string>
#include <iostream>
#include <vector>

#include "data-struct/Bipartition.hpp"

// class Bipartition; 
typedef unsigned int nat; 

#ifdef BYTE_ALIGNMENT
#undef BYTE_ALIGNMENT
#endif
#define BYTE_ALIGNMENT 32

namespace  Parser
{
#include "parser/axml-parser.hpp"
}
#undef BYTE_ALIGNMENT

class PhylipParser
{
public: 
  PhylipParser(std::string _alnFile, std::string _modelFile, bool _haveModelFile); 
  ~PhylipParser();
  void parse(); 
  void writeToFile(std::string fileName) ; 

private: 			// METHODS
  void writeWeights(std::ofstream &out); 

  template<typename T>
  void myWrite(std::ostream& out, T* ptr, size_t length)
  {
    out.write((char*) ptr, sizeof(T) * length) ; 
  }

  void getyspace();
  int iterated_bitcount(unsigned int n);
  void defaultInit();
  void parseHeader(FILE *INFILE);
  Parser::boolean setupTree ();
  void checkTaxonName(char *buffer, int len);
  void parsePartitions(); 
  void parseSinglePartition(std::string dataType); 
  void uppercase (int *chptr);
  Parser::boolean getdata(FILE *INFILE);
  unsigned char buildStates(int secModel, unsigned char v1, unsigned char v2);
  void adaptRdataToSecondary();
  void sitesort();
  void sitecombcrunch ();
  Parser::boolean makeweights();
  Parser::boolean makevalues();
  Parser::hashNumberType  hashString(char *p, Parser::hashNumberType tableSize);
  void addword(char *s, Parser::stringHashtable *h, int nodeNumber);
  Parser::stringHashtable* initStringHashTable(Parser::hashNumberType n);
  void getinput();

  int getStates(int dataType); 
  int getUndetermined(int dataType); 
  const Parser::partitionLengths* getPartitionLengths(Parser::pInfo *p); 
  const unsigned int* getBitVector(int dataType); 

  // parsimony stuff 
  void compressDNA(int *informative); 
  Parser::boolean isInformative(int dataType, int site); 
  void determineUninformativeSites( int *informative); 
  Bipartition allocateParsimonyDataStructures(); 

private: 			// ATTRIBUTES
  std::string alnFile; 		
  std::string modelFile; 

  Parser::rawdata      *rdta;
  Parser::cruncheddata *cdta;
  Parser::tree         *tr;
  Parser::analdef      *adef;

  bool _compress; 		// an option 
  nat _numTax; 
  nat _numSites; 

  unsigned int myMask32[32]; 
  std::vector<std::string> protModels; 

  const char inverseMeaningBINARY[4] ; 
  const char inverseMeaningDNA[16]; 
  const char inverseMeaningPROT[23]; 
  const char inverseMeaningGeneric32[33]; 
  const char inverseMeaningGeneric64[33]; 
  const unsigned int bitVectorIdentity[256]; 
  const unsigned int bitVectorAA[23]; 
  const unsigned int bitVectorSecondary[256]; 
  const unsigned int bitVector32[33]; 
  Parser::partitionLengths pLengths[8]; 
  
  bool _haveModelFile; 

  Parser::nodeptr baseAddr; 
  
  std::vector<std::vector<uint8_t> > _ySpace; // 

  Bipartition _infoness; 
}; 


#endif


