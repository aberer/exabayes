#include "parser/PhylipParser.hpp"

/**
   @brief disclaimer: this is a relatively quick and very dirty
   approach to integrating the parser into the main code
 */ 

#include <assert.h>

#include "../system/time.hpp"

#include "data-struct/Bipartition.hpp"

#include "../common.h"
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>


// #define OLD_ALN_LAYOUT

// #define DEBUG_MSG


using namespace Parser; 

#include "parser/parserDefines.hpp"

#define INTS_PER_VECTOR 8

typedef unsigned int  nat; 
extern const char *protModels[NUM_PROT_MODELS];


static void myExit(int code)
{
  exit(code); 
}


PhylipParser::PhylipParser(std::string _alnFile, std::string _modelFile, bool haveModelFile	)
  : alnFile(_alnFile)
  , modelFile(_modelFile)
  , _compress(true)		// TODO 
  , protModels {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM", "LG", "MTART", "MTZOA", "PMB", "HIVB", "HIVW", "JTTDCMUT", "FLU", "AUTO", "LG4", "GTR"}
  , inverseMeaningBINARY {'_', '0', '1', '-'}
  , inverseMeaningDNA {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'}
  , inverseMeaningPROT {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', '-'}
  , inverseMeaningGeneric32 {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', '-'}
  , inverseMeaningGeneric64 {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', '-'}
  , bitVectorIdentity  {0 ,1 ,2 ,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22 ,23 ,24 ,25 ,26 , 27 ,28 ,29 ,30 ,31 ,32 ,33 ,34 ,35 ,36 ,37 ,38 ,39 ,40 ,41 ,42 ,43 ,44 ,45 ,46 ,47 ,48 ,49 ,50 ,51 , 52 ,53 ,54 ,55 ,56 ,57 ,58 ,59 ,60 ,61 ,62 ,63 ,64 ,65 ,66 ,67 ,68 ,69 ,70 ,71 ,72 ,73 ,74 ,75 ,76 , 77 ,78 ,79 ,80 ,81 ,82 ,83 ,84 ,85 ,86 ,87 ,88 ,89 ,90 ,91 ,92 ,93 ,94 ,95 ,96 ,97 ,98 ,99 ,100 ,101 , 102 ,103 ,104 ,105 ,106 ,107 ,108 ,109 ,110 ,111 ,112 ,113 ,114 ,115 ,116 ,117 ,118 ,119 ,120 ,121 ,122 , 123 ,124 ,125 ,126 ,127 ,128 ,129 ,130 ,131 ,132 ,133 ,134 ,135 ,136 ,137 ,138 ,139 ,140 ,141 ,142 ,143 , 144 ,145 ,146 ,147 ,148 ,149 ,150 ,151 ,152 ,153 ,154 ,155 ,156 ,157 ,158 ,159 ,160 ,161 ,162 ,163 ,164 , 165 ,166 ,167 ,168 ,169 ,170 ,171 ,172 ,173 ,174 ,175 ,176 ,177 ,178 ,179 ,180 ,181 ,182 ,183 ,184 ,185 , 186 ,187 ,188 ,189 ,190 ,191 ,192 ,193 ,194 ,195 ,196 ,197 ,198 ,199 ,200 ,201 ,202 ,203 ,204 ,205 ,206 , 207 ,208 ,209 ,210 ,211 ,212 ,213 ,214 ,215 ,216 ,217 ,218 ,219 ,220 ,221 ,222 ,223 ,224 ,225 ,226 ,227 , 228 ,229 ,230 ,231 ,232 ,233 ,234 ,235 ,236 ,237 ,238 ,239 ,240 ,241 ,242 ,243 ,244 ,245 ,246 ,247 ,248 , 249 ,250 ,251 ,252 ,253 ,254 ,255}
  , bitVectorAA {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 12 /* N | D */, 96 /*Q | E*/, 1048575 /* - */}
  ,  bitVectorSecondary {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240, 0, 17, 34, 51, 68, 85, 102, 119, 136, 153, 170, 187, 204, 221, 238, 255, 0, 256, 512, 768, 1024, 1280, 1536, 1792, 2048, 2304, 2560, 2816, 3072, 3328, 3584, 3840, 0, 257, 514, 771, 1028, 1285, 1542, 1799, 2056, 2313, 2570, 2827, 3084, 3341, 3598, 3855, 0, 272, 544, 816, 1088, 1360, 1632, 1904, 2176, 2448, 2720, 2992, 3264, 3536, 3808, 4080, 0, 273, 546, 819, 1092, 1365, 1638, 1911, 2184, 2457, 2730, 3003, 3276, 3549, 3822, 4095, 0, 4096, 8192, 12288, 16384, 20480, 24576, 28672, 32768, 36864, 40960, 45056, 49152, 53248, 57344, 61440, 0, 4097, 8194, 12291, 16388, 20485, 24582, 28679, 32776, 36873, 40970, 45067, 49164, 53261, 57358, 61455, 0, 4112, 8224, 12336, 16448, 20560, 24672, 28784, 32896, 37008, 41120, 45232, 49344, 53456, 57568, 61680, 0, 4113, 8226, 12339, 16452, 20565, 24678, 28791, 32904, 37017, 41130, 45243, 49356, 53469, 57582, 61695, 0, 4352, 8704, 13056, 17408, 21760, 26112, 30464, 34816, 39168, 43520, 47872, 52224, 56576, 60928, 65280, 0, 4353, 8706, 13059, 17412, 21765, 26118, 30471, 34824, 39177, 43530, 47883, 52236, 56589, 60942, 65295, 0, 4368, 8736, 13104, 17472, 21840, 26208, 30576, 34944, 39312, 43680, 48048, 52416, 56784, 61152, 65520, 0, 4369, 8738, 13107, 17476, 21845, 26214, 30583, 34952, 39321, 43690, 48059, 52428, 56797, 61166, 65535}
  ,  bitVector32 {1,     2,    4,    8,   16,   32,    64,   128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648u, 4294967295u}
  , pLengths  {
    {4,   4,   2,  4,  2, 1, 2,  8, 2, 2, FALSE, FALSE, 3, inverseMeaningBINARY, 2, FALSE, bitVectorIdentity},    /* BINARY */
    {16,  16,  4, 16, 16, 6, 4, 64, 6, 4, FALSE, FALSE, 15, inverseMeaningDNA, 4, FALSE, bitVectorIdentity},    /* DNA */
    {400, 400, 20, 400, 400, 190, 20, 460, 190, 20, FALSE, FALSE, 22, inverseMeaningPROT, 20, TRUE, bitVectorAA},    /* AA */
    {256, 256, 16, 256, 256, 120, 16, 4096, 120, 16, FALSE, FALSE, 255, (char*)NULL, 16, TRUE, bitVectorSecondary},    /* SECONDARY_DATA */
    {36, 36,  6, 36, 36, 15, 6, 384, 15, 6, FALSE, FALSE, 63, (char*)NULL, 6, TRUE, bitVectorIdentity},    /* SECONDARY_DATA_6 */
    {49,   49,    7,   49, 49,  21, 7, 896, 21, 7, FALSE, FALSE, 127, (char*)NULL, 7, TRUE, bitVectorIdentity},    /* SECONDARY_DATA_7 */
    {1024, 1024, 32, 1024, 1024, 496, 32, 1056, 496, 32, FALSE, FALSE, 32, inverseMeaningGeneric32, 32, TRUE, bitVector32},    /* 32 states */
    {4096, 4096, 64, 4096, 4096, 2016, 64, 4160, 64, 2016, FALSE, FALSE, 64, (char*)NULL, 64, TRUE, (unsigned int*)NULL}    /* 64 states */
  }
  , _haveModelFile(haveModelFile)
{
  // std::cout << "called parser with " << alnFile << "\t" << modelFile << " and " << ( haveModelFile ? "have model" : "NO model"  ) <<  std::endl; 
  
  myMask32[0] = 1 ; 
  for(nat i = 1; i < 32; ++i)
    myMask32[i] = myMask32[i-1] << 1 ; 
}


PhylipParser::~PhylipParser()
{
  // if(rdta && rdta->y0)
  free(rdta->y0); 
  free(baseAddr); 
  free(tr->yVector); 
  free(tr->ti); 

  free(tr->nameHash->table); 
  free(tr->nameHash);


  free(rdta->wgt); 
  free(cdta->alias); 
  free(cdta->aliaswgt); 
  free(tr->model); 
  free(tr->initialDataVector); 
  free(tr->extendedDataVector); 

    for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      auto *partition = tr->partitionData + i; 
      free(partition->parsVect); 
      free(partition->partitionName); 
    }

  for(int i = 1; i < tr->mxtips + 1; ++i)  
    free(tr->nameList[i]); 
  free(tr->nameList); 

  free(tr->initialPartitionData); 
  free(tr); 
  free(cdta); 
  free(rdta); 
  free(adef); 
}



int PhylipParser::getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);
  return pLengths[dataType].states;
}

int PhylipParser::getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);
  return pLengths[dataType].undetermined;
}


static boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


const partitionLengths* PhylipParser::getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return (&pLengths[dataType]); 
}


static void skipWhites(char **ch)
{
  while(**ch == ' ' || **ch == '\t')
    *ch = *ch + 1;
}


static boolean lineContainsOnlyWhiteChars(char *line)
{
  int i, n = strlen(line);

  if(n == 0)
    return TRUE;

  for(i = 0; i < n; i++)
    {
      if(!whitechar(line[i]))
	return FALSE;
    }
  return TRUE;
}


const unsigned int* PhylipParser::getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


static int myGetline(char **lineptr, int *n, FILE *stream)
{
  char *line, *p;
  int size, copy, len;
  int chunkSize = 256 * sizeof(char);

  if (*lineptr == NULL || *n < 2) 
    {
      line = (char *)realloc(*lineptr, chunkSize);
      if (line == NULL)
	return -1;
      *lineptr = line;
      *n = chunkSize;
    }

  line = *lineptr;
  size = *n;
  
  copy = size;
  p = line;
   
  while(1)
    {
      while (--copy > 0)
	{
	  int c = getc(stream);
	  if (c == EOF)
	    goto lose;
	  else
	    {
	      *p++ = c;
	      if(c == '\n' || c == '\r')	
		goto win;
	    }
	}

      /* Need to enlarge the line buffer.  */
      len = p - line;
      size *= 2;
      line = (char*) realloc (line, size);
      if (line == NULL)
	goto lose;
      *lineptr = line;
      *n = size;
      p = line + len;
      copy = size - len;
    }
   
 lose:
  if (p == *lineptr)
    return -1;
  /* Return a partial line since we got an error in the middle.  */
 win:
  *p = '\0';
  return p - *lineptr;
}

  

static void analyzeIdentifier(char **ch, int modelNumber, tree *tr)
{
  char ident[2048] = "";
  char model[128] = "";  
  int  j;
  int containsComma = 0;
  
  int n = 0; 
  while(**ch != '=')
    {
      if(**ch != ' ' && **ch != '\t')
	{
	  ident[n] = **ch;      
	  n++;
	}
      *ch = *ch + 1;
    }

  for(int i = 0; i < n; i++)
    if(ident[i] == ',') 
      containsComma = 1;

  if(!containsComma)
    {
      printf("Error, model file must have format: DNA or AA model, then a comma, and then the partition name\n");
      myExit(-1);
    }

  int i = 0;
  while(ident[i] != ',')
    {
      model[i] = ident[i];
      i++;
    }      

  /* AA */
  if( strcasecmp(model, "PROT") == 0 )
    {
      tr->initialPartitionData[modelNumber].protModels = i;		  
      tr->initialPartitionData[modelNumber].protFreqs  = 0;
      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;
      tr->initialPartitionData[modelNumber].dataType   = AA_DATA;
    }
  else if(strcasecmp(model, "DNA") == 0)
    {	     	      
      tr->initialPartitionData[modelNumber].protModels = -1;		  
      tr->initialPartitionData[modelNumber].protFreqs  = -1;
      tr->initialPartitionData[modelNumber].dataType   = DNA_DATA;
      tr->initialPartitionData[modelNumber].optimizeBaseFrequencies = FALSE;
    }
  else if (strcasecmp(model, "BIN") == 0)
    {	     	      
      std::cout << "model >BIN< is not supported yet." << std::endl; 
      myExit(-1); 
	
      tr->initialPartitionData[modelNumber].protModels = -1;		  
      tr->initialPartitionData[modelNumber].protFreqs  = -1;
      tr->initialPartitionData[modelNumber].dataType   = BINARY_DATA;
    }
  else if(strcasecmp(model, "MULTI") == 0)
    {	     	
      std::cout << "model >MULTI< is not supported yet." << std::endl; 
      myExit(-1); 
      
      tr->initialPartitionData[modelNumber].protModels = -1;		  
      tr->initialPartitionData[modelNumber].protFreqs  = -1;
      tr->initialPartitionData[modelNumber].dataType   = GENERIC_32;
    }
  else if(strcasecmp(model, "CODON") == 0)
    {	     	      
      std::cout << "model >CODON< is not supported yet." << std::endl; 
      myExit(-1); 

      tr->initialPartitionData[modelNumber].protModels = -1;		  
      tr->initialPartitionData[modelNumber].protFreqs  = -1;
      tr->initialPartitionData[modelNumber].dataType   = GENERIC_64;
    }

  i = 0; 
  while(ident[i++] != ',');      
  tr->initialPartitionData[modelNumber].partitionName = (char*)malloc((n - i + 1) * sizeof(char));          
  j = 0;
  while(i < n)	
    tr->initialPartitionData[modelNumber].partitionName[j++] =  ident[i++];
  tr->initialPartitionData[modelNumber].partitionName[j] = '\0';                      
}


static int isNum(char c)
{
  
  return (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
	  c == '5' || c == '6' || c == '7' || c == '8' || c == '9');
 }

static void setModel(int model, int position, int *a)
{
  if(a[position] == -1)
    a[position] = model;
  else
    {
      printf("ERROR trying to assign model %d to position %d \n", model, position);
      printf("while already model %d has been assigned to this position\n", a[position]);
      myExit(-1);
    }      
}



void PhylipParser::parseSinglePartition(std::string dataType)
{
  tr->NumberOfModels = 1; 
  tr->initialPartitionData = (pInfo*)calloc(1,sizeof(pInfo));
  
  auto partition = tr->initialPartitionData; 
  
  partition->protModels = adef->proteinMatrix;
  partition->protFreqs  = adef->protEmpiricalFreqs;

  partition->lower = 0; 
  partition->upper = tr->originalCrunchedLength; 
  partition->width = tr->originalCrunchedLength; 
  
  if(dataType.compare("PROT") == 0)
    {
      partition->dataType = AA_DATA; 
      partition->states =  20; 
    }
  else if(dataType.compare("DNA") == 0)
    {
      partition->dataType = DNA_DATA; 
      partition->states = 4; 

    }
  else 
    {
      std::cout << "parameter for -m must be either PROT or DNA" << std::endl; 
      assert(0); 
    }

  partition->partitionName = strdup("NoName"); 

}

void PhylipParser::parsePartitions()
{
  FILE *f; 
  nat numberOfModels = 0; 
  int nbytes = 0;
  char *ch;
  char *cc = (char*)calloc(1,sizeof(char));
  char **p_names;
  int n, l;
  int lower, upper, modulo;
  char buf[256];
  int **partitions;
  int pairsCount;
  nat as;
  // int k; 

  f = fopen(modelFile.c_str(), "rb");   

  while(myGetline(&cc, &nbytes, f) > -1)
    {     
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  numberOfModels++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }     

  rewind(f);
      
  p_names = (char **)malloc(sizeof(char *) * numberOfModels);
  partitions = (int **)malloc(sizeof(int *) * numberOfModels);
  
  tr->initialPartitionData = (pInfo*)calloc(numberOfModels, sizeof(pInfo)  );

  for(nat i = 0; i < numberOfModels; i++) 
    {     
      tr->initialPartitionData[i].protModels = adef->proteinMatrix;
      tr->initialPartitionData[i].protFreqs  = adef->protEmpiricalFreqs;
      tr->initialPartitionData[i].dataType   = -1;
    }

  for(nat i = 0; i < numberOfModels; i++)    
    partitions[i] = (int *)NULL;

  nat i = 0; 
  while(myGetline(&cc, &nbytes, f) > -1)
    {          
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  n = strlen(cc);	 
	  p_names[i] = (char *)malloc(sizeof(char) * (n + 1));
	  strcpy(&(p_names[i][0]), cc);
	  i++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }         

  for(i = 0; i < numberOfModels; i++)
    {           
      ch = p_names[i];     
      pairsCount = 0;
      skipWhites(&ch);
      
      if(*ch == '=')
	{
	  printf("Identifier missing prior to '=' in %s\n", p_names[i]);
	  myExit(-1);
	}
      
      analyzeIdentifier(&ch, i, tr);
      ch++;
            
    numberPairs:
      pairsCount++;
      partitions[i] = (int *)realloc((void *)partitions[i], (1 + 3 * pairsCount) * sizeof(int));
      partitions[i][0] = pairsCount;
      partitions[i][3 + 3 * (pairsCount - 1)] = -1; 	
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  myExit(-1);
	}   
      
      l = 0;
      while(isNum(*ch))		 
	{
	  /*printf("%c", *ch);*/
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      lower = atoi(buf);
      partitions[i][1 + 3 * (pairsCount - 1)] = lower;   
      
      skipWhites(&ch);
      
      /* NEW */
      
      if((*ch != '-') && (*ch != ','))
	{
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {
	      upper = lower;
	      goto SINGLE_NUMBER;
	    }
	  else
	    {
	      printf("'-' or ',' expected in %s\n", p_names[i]);
	      myExit(-1);
	    }
	}	 
      
      if(*ch == ',')
	{	     
	  upper = lower;
	  goto SINGLE_NUMBER;
	}
      
      /* END NEW */
      
      ch++;   
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  myExit(-1);
	}    
      
      l = 0;
      while(isNum(*ch))
	{    
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      upper = atoi(buf);     
    SINGLE_NUMBER:
      partitions[i][2 + 3 * (pairsCount - 1)] = upper;        	  
      
      if(upper < lower)
	{
	  printf("Upper bound %d smaller than lower bound %d for this partition: %s\n", upper, lower,  p_names[i]);
	  myExit(-1);
	}
      
      skipWhites(&ch);
      
      if(*ch == '\0' || *ch == '\n' || *ch == '\r') /* PC-LINEBREAK*/
	{    
	  goto parsed;
	}
      
      if(*ch == ',')
	{	 
	  ch++;
	  goto numberPairs;
	}
      
      if(*ch == '\\')
	{
	  ch++;
	  skipWhites(&ch);
	  
	  if(!isNum(*ch))
	    {
	      printf("%c Number expected in %s\n", *ch, p_names[i]);
	      myExit(-1);
	    }     
	  
	  l = 0;
	  while(isNum(*ch))
	    {
	      buf[l] = *ch;
	      ch++;	
	      l++;
	    }
	  buf[l] = '\0';
	  modulo = atoi(buf);      
	  partitions[i][3 + 3 * (pairsCount - 1)] = modulo; 	
	  
	  skipWhites(&ch);
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {	     
	      goto parsed;
	    }
	  if(*ch == ',')
	    {	       
	      ch++;
	      goto numberPairs;
	    }
	}  
      
      assert(0);
       
    parsed:
      ; 
      // i = i; 			// ??? 
    }
  
  fclose(f);
 
  /*********************************************************************************************************************/ 

  for(nat i = 0; i <= _numSites; i++)
    tr->model[i] = -1;
  
  for(nat i = 0; i < numberOfModels; i++)
    {   
      as = partitions[i][0];     
      
      for(nat j = 0; j < as; j++)
	{
	  lower = partitions[i][1 + j * 3];
	  upper = partitions[i][2 + j * 3]; 
	  modulo = partitions[i][3 + j * 3];	
	 
	  if(modulo == -1)
	    {
	      for(int k = lower; k <= upper; k++)
		setModel(i, k, tr->model);
	    }
	  else
	    {
	      for(int k = lower; k <= upper; k += modulo)
		{
		  if(k <= int(_numSites))
		    setModel(i, k, tr->model);	      
		}
	    }
	}        
    }


  for(i = 1; i < _numSites + 1; i++)
    {
      
      if(tr->model[i] == -1)
	{
	  printf("ERROR: Alignment Position %d has not been assigned any model\n", i);
	  myExit(-1);
	}      
    }  

  for(i = 0; i < numberOfModels; i++)
    {
      free(partitions[i]);
      free(p_names[i]);
    }
  
  free(partitions);
  free(p_names);    
    
  tr->NumberOfModels = numberOfModels;     
  
  if(adef->perGeneBranchLengths)
    tr->numBranches = tr->NumberOfModels;
}

void PhylipParser::getyspace ()
{
  size_t size = 4 * ((size_t)(_numSites / 4 + 1));
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc((_numTax + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *) malloc(((size_t)(_numTax + 1)) * size * sizeof(unsigned char));

  /*
    printf("Raw alignment data Assigning %Zu bytes\n", ((size_t)(_numTax + 1)) * size * sizeof(unsigned char));
  */

  assert(y0);   

  rdta->y0 = y0;

  for (nat i = 0; i <= _numTax; i++)
    {
      rdta->y[i] = y0;
      y0 += size;
    }
  return;
}


void PhylipParser::defaultInit()
{  
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25;
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE;
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;  
  adef->useInvariant           = FALSE;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE;
  adef->allInOne               = FALSE;
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->externalAAMatrix       = (double*)NULL;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->thoroughInsertion      = FALSE;
  adef->compressPatterns       = TRUE; 
  adef->readTaxaOnly           = FALSE;
  adef->meshSearch             = 0;
  adef->useCheckpoint          = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->catOnly = FALSE;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;

  if(_compress)
    adef->compressPatterns = TRUE; 
  
  adef->useMultipleModel = TRUE; 
}


void PhylipParser::parseHeader (FILE *INFILE)
{
   if (fscanf(INFILE, "%u %u", & _numTax, & _numSites) != 2)
    {
      
      printf("\n Error: problem reading number of species and sites\n\n");
      myExit(-1);
    }
  
  if (_numTax < 4) // _numTax
    {
      printf("\n Error: too few species\n\n");
      myExit(-1);
    }

  if (_numSites < 1) // _numSites
    {
      printf("\n Error: too few sites\n\n");
      myExit(-1);
    }

  return;
}


boolean PhylipParser::setupTree ()
{
  nodeptr  
    p0;
  
  int
    tips,
    inter; 

  if(!adef->readTaxaOnly)
    {
      tr->patternPosition = (int*)NULL;
      tr->columnPosition = (int*)NULL;
    }

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  if(!adef->readTaxaOnly)
    tr->yVector      = (unsigned char **)  malloc((tr->mxtips + 1) * sizeof(unsigned char *));
    
  /*TODO, must that be so long ?*/

  if(!adef->readTaxaOnly)
    {
      tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));
    }

  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("\n Error: unable to obtain sufficient tree memory\n\n");
      return  FALSE;
    }

  baseAddr = p0; 
  
  tr->vLength = 0;
  tr->h = (hashtable*)NULL;

  return TRUE;
}


void PhylipParser::checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':
	case '\t':
	case '\n':
	case '\r':
	case ' ':
	case ':':
	case ',':
	case '(':
	case ')':
	case ';':
	case '[':
	case ']':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}

      if(!valid)
	{
	  printf("\n Error: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n\n", buffer, i, buffer[i]);
	  printf(" Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\"\n");
	  printf(" MyExiting\n");
	  myExit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}

void PhylipParser::uppercase (int *chptr)
{
  int  ch;

  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}

boolean PhylipParser::getdata(FILE *INFILE)
{
  int   
    i, 
    j, 
    basesread, 
    basesnew, 
    ch, my_i, meaning,
    len,
    meaningAA[256], 
    meaningDNA[256], 
    meaningBINARY[256],
    meaningGeneric32[256],
    meaningGeneric64[256];
  
  boolean  
    allread, 
    firstpass;
  
  char 
    buffer[nmlngth + 2];
  
  unsigned char
    genericChars32[32] = {'0', '1', '2', '3', '4', '5', '6', '7', 
			  '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
			  'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
			  'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};  
  unsigned long 
    total = 0,
    gaps  = 0;

  for (i = 0; i < 256; i++)
    {      
      meaningAA[i]          = -1;
      meaningDNA[i]         = -1;
      meaningBINARY[i]      = -1;
      meaningGeneric32[i]   = -1;
      meaningGeneric64[i]   = -1;
    }

  /* generic 32 data */

  for(i = 0; i < 32; i++)
    meaningGeneric32[genericChars32[i]] = i;
  meaningGeneric32[int('-')] = getUndetermined(GENERIC_32);
  meaningGeneric32[int('?')] = getUndetermined(GENERIC_32);

  /* AA data */

  meaningAA[int('A')] =  0;  /* alanine */
  meaningAA[int('R')] =  1;  /* arginine */
  meaningAA[int('N')] =  2;  /*  asparagine*/
  meaningAA[int('D')] =  3;  /* aspartic */
  meaningAA[int('C')] =  4;  /* cysteine */
  meaningAA[int('Q')] =  5;  /* glutamine */
  meaningAA[int('E')] =  6;  /* glutamic */
  meaningAA[int('G')] =  7;  /* glycine */
  meaningAA[int('H')] =  8;  /* histidine */
  meaningAA[int('I')] =  9;  /* isoleucine */
  meaningAA[int('L')] =  10; /* leucine */
  meaningAA[int('K')] =  11; /* lysine */
  meaningAA[int('M')] =  12; /* methionine */
  meaningAA[int('F')] =  13; /* phenylalanine */
  meaningAA[int('P')] =  14; /* proline */
  meaningAA[int('S')] =  15; /* serine */
  meaningAA[int('T')] =  16; /* threonine */
  meaningAA[int('W')] =  17; /* tryptophan */
  meaningAA[int('Y')] =  18; /* tyrosine */
  meaningAA[int('V')] =  19; /* valine */
  meaningAA[int('B')] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA[int('Z')] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA[int('X')] = 
    meaningAA[int('?')] = 
    meaningAA[int('*')] = 
    meaningAA[int('-')] = 
    getUndetermined(AA_DATA);

  /* DNA data */

  meaningDNA[int('A')] =  1;
  meaningDNA[int('B')] = 14;
  meaningDNA[int('C')] =  2;
  meaningDNA[int('D')] = 13;
  meaningDNA[int('G')] =  4;
  meaningDNA[int('H')] = 11;
  meaningDNA[int('K')] = 12;
  meaningDNA[int('M')] =  3;  
  meaningDNA[int('R')] =  5;
  meaningDNA[int('S')] =  6;
  meaningDNA[int('T')] =  8;
  meaningDNA[int('U')] =  8;
  meaningDNA[int('V')] =  7;
  meaningDNA[int('W')] =  9; 
  meaningDNA[int('Y')] = 10;

  meaningDNA[int('N')] = meaningDNA[int('O')] = meaningDNA[int('X')] = meaningDNA[int('-')] = meaningDNA[int('?')] = getUndetermined(DNA_DATA);

  /* BINARY DATA */

  meaningBINARY[int('0')] = 1;
  meaningBINARY[int('1')] = 2;
  
  meaningBINARY[int('-')] = 
    meaningBINARY[int('?')] = 
    getUndetermined(BINARY_DATA);


  /*******************************************************************/

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
    {
      for (i = 1; i <= tr->mxtips; i++)
	{
	  if (firstpass)
	    {
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);

	      my_i = 0;

	      do
		{
		  buffer[my_i] = ch;
		  ch = getc(INFILE);
		  my_i++;
		  if(my_i >= nmlngth)
		    {

		      printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
		      printf("axml.h, current setting %d\n", nmlngth);

		      myExit(-1);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      
	      ungetc(ch, INFILE);

	      buffer[my_i] = '\0';
	      len = strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * len);
	      strcpy(tr->nameList[i], buffer);
	    }

	  j = basesread;

	  while ((j < int(_numSites)) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {
	      uppercase(& ch);

	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case BINARY_DATA:
		  meaning = meaningBINARY[ch];
		  break;
		case DNA_DATA:
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7:
		  /*
		    still dealing with DNA/RNA here, hence just act if as they where DNA characters
		    corresponding column merging for sec struct models will take place later
		  */
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		case GENERIC_32:
		  meaning = meaningGeneric32[ch];
		  break;
		case GENERIC_64:
		  meaning = meaningGeneric64[ch];
		  break;
		default:
		  assert(0);
		}

	      if (meaning != -1)
		{
		  j++;
		  rdta->y[i][j] = ch;		 
		}
	      else
		{
		  if(!whitechar(ch))
		    {
		      printf("\n Error: bad base (%c) at site %d of sequence %d\n\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
	    }

	  if (ch == EOF)
	    {
	      printf("\n Error: end-of-file at site %d of sequence %d\n\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread))
	    i--;
	  else
	    {
	      if (i == 1)
		basesnew = j;
	      else
		if (j != basesnew)
		  {
		    printf("\n Error: sequences out of alignment\n");
		    printf("%d (instead of %d) residues read in sequence %d %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= int(_numSites));
    }

  for(j = 1; j <= tr->mxtips; j++)
    for(i = 1; i <= int(_numSites); i++)
      {
	assert(tr->dataVector[i] != -1);

	switch(tr->dataVector[i])
	  {
	  case BINARY_DATA:
	    meaning = meaningBINARY[rdta->y[j][i]];
	    if(meaning == getUndetermined(BINARY_DATA))
	      gaps++;
	    break;

	  case SECONDARY_DATA:
	  case SECONDARY_DATA_6:
	  case SECONDARY_DATA_7:
	    assert(tr->secondaryStructurePairs[i - 1] != -1);
	    assert(i - 1 == tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]]);
	    /*
	      don't worry too much about undetermined column count here for sec-struct, just count
	      DNA/RNA gaps here and worry about the rest later-on, falling through to DNA again :-)
	    */
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == getUndetermined(DNA_DATA))
	      gaps++;
	    break;

	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == getUndetermined(AA_DATA))
	      gaps++;
	    break;

	  case GENERIC_32:
	    meaning = meaningGeneric32[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_32))
	      gaps++;
	    break;

	  case GENERIC_64:
	    meaning = meaningGeneric64[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_64))
	      gaps++;
	    break;
	  default:
	    assert(0);
	  }

	total++;
	rdta->y[j][i] = meaning;
      }

  adef->gapyness = (double)gaps / (double)total;

  return  TRUE;
}


unsigned char PhylipParser::buildStates(int secModel, unsigned char v1, unsigned char v2)
{
  unsigned char newState = 0;

  switch(secModel)
    {
    case SECONDARY_DATA:
      newState = v1;
      newState = newState << 4;
      newState = newState | v2;
      break;
    case SECONDARY_DATA_6:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[6] = {1, 2, 4, 8, 16, 32};

	unsigned char
	  intermediateBinaryStates[6];

	int length = 6;

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA[int('A')] =  1;
	meaningDNA[int('B')] = 14;
	meaningDNA[int('C')] =  2;
	meaningDNA[int('D')] = 13;
	meaningDNA[int('G')] =  4;
	meaningDNA[int('H')] = 11;
	meaningDNA[int('K')] = 12;
	meaningDNA[int('M')] =  3;
	meaningDNA[int('N')] = 15;
	meaningDNA[int('O')] = 15;
	meaningDNA[int('R')] =  5;
	meaningDNA[int('S')] =  6;
	meaningDNA[int('T')] =  8;
	meaningDNA[int('U')] =  8;
	meaningDNA[int('V')] =  7;
	meaningDNA[int('W')] =  9;
	meaningDNA[int('X')] = 15;
	meaningDNA[int('Y')] = 10;
	meaningDNA[int('-')] = 15;
	meaningDNA[int('?')] = 15;

	for(i = 0; i < length; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    newState = n1;
	    newState = newState << 4;
	    newState = newState | n2;

	    intermediateBinaryStates[i] = newState;
	  }

	newState = v1;
	newState = newState << 4;
	newState = newState | v2;

	for(i = 0; i < length; i++)
	  {
	    if(newState == intermediateBinaryStates[i])
	      break;
	  }
	if(i < length)
	  newState = finalBinaryStates[i];
	else
	  {
	    newState = 0;
	    for(i = 0; i < length; i++)
	      {
		if(v1 & meaningDNA[allowedStates[i][0]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    newState |= finalBinaryStates[i];
		  }
		if(v2 & meaningDNA[allowedStates[i][1]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    newState |= finalBinaryStates[i];
		  }
	      }
	  }	
      }
      break;
    case SECONDARY_DATA_7:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[7] = {1, 2, 4, 8, 16, 32, 64};

	unsigned char
	  intermediateBinaryStates[7];

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA[int('A')] =  1;
	meaningDNA[int('B')] = 14;
	meaningDNA[int('C')] =  2;
	meaningDNA[int('D')] = 13;
	meaningDNA[int('G')] =  4;
	meaningDNA[int('H')] = 11;
	meaningDNA[int('K')] = 12;
	meaningDNA[int('M')] =  3;
	meaningDNA[int('N')] = 15;
	meaningDNA[int('O')] = 15;
	meaningDNA[int('R')] =  5;
	meaningDNA[int('S')] =  6;
	meaningDNA[int('T')] =  8;
	meaningDNA[int('U')] =  8;
	meaningDNA[int('V')] =  7;
	meaningDNA[int('W')] =  9;
	meaningDNA[int('X')] = 15;
	meaningDNA[int('Y')] = 10;
	meaningDNA[int('-')] = 15;
	meaningDNA[int('?')] = 15;
	

	for(i = 0; i < 6; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    newState = n1;
	    newState = newState << 4;
	    newState = newState | n2;

	    intermediateBinaryStates[i] = newState;
	  }

	newState = v1;
	newState = newState << 4;
	newState = newState | v2;

	for(i = 0; i < 6; i++)
	  {
	    /* exact match */
	    if(newState == intermediateBinaryStates[i])
	      break;
	  }
	if(i < 6)
	  newState = finalBinaryStates[i];
	else
	  {
	    /* distinguish between exact mismatches and partial mismatches */

	    for(i = 0; i < 6; i++)
	      if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		break;
	    if(i < 6)
	      {
		/* printf("partial mismatch\n"); */

		newState = 0;
		for(i = 0; i < 6; i++)
		  {
		    if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		      {
			/*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
			newState |= finalBinaryStates[i];
		      }
		    else
		      newState |=  finalBinaryStates[6];
		  }
	      }
	    else
	      newState = finalBinaryStates[6];
	  }	
      }
      break;
    default:
      assert(0);
    }

  return newState;
}


void PhylipParser::adaptRdataToSecondary()
{
  int *alias = (int*)calloc(_numSites, sizeof(int));
  int  realPosition = 0;  

  for(nat i = 0; i < _numSites; i++)
    alias[i] = -1;

  for(nat i = 0; i < _numSites; i++)
    {
      int partner = tr->secondaryStructurePairs[i];
      if(partner != -1)
	{
	  assert(tr->dataVector[i+1] == SECONDARY_DATA || tr->dataVector[i+1] == SECONDARY_DATA_6 || tr->dataVector[i+1] == SECONDARY_DATA_7);

	  if(int(i) < partner)
	    {
	      for(nat j = 1; j <= _numTax; j++)
		{
		  unsigned char v1 = rdta->y[j][i+1];
		  unsigned char v2 = rdta->y[j][partner+1];

		  assert(int(i+1) < partner+1);

		  rdta->y[j][i+1] = buildStates(tr->dataVector[i+1], v1, v2);
		}
	      alias[realPosition] = i;
	      realPosition++;
	    }
	}
      else
	{
	  alias[realPosition] = i;
	  realPosition++;
	}
    }

  assert(int(_numSites - realPosition) == tr->numberOfSecondaryColumns / 2);

  _numSites = realPosition;

  for(nat i = 0; i < _numSites; i++)
    {
      assert(alias[i] != -1);
      tr->model[i+1]    = tr->model[alias[i]+1];
      tr->dataVector[i+1] = tr->dataVector[alias[i]+1];
      rdta->wgt[i+1] =  rdta->wgt[alias[i]+1];

      for(nat j = 1; j <= _numTax; j++)
	rdta->y[j][i+1] = rdta->y[j][alias[i]+1];
    }

  free(alias);
}


// #define OLD_SORT 

void PhylipParser::sitesort()
{
  unsigned char  **data;
  int  
    *index, 
    *category = (int*)NULL;
 
#ifdef OLD_SORT
  int n = _numSites, 
    nsp = _numTax; 
  int  gap, i, j, jj, jg, k;
  boolean  flip, tied;
#endif

  if(adef->useSecondaryStructure)
    {
      assert(tr->NumberOfModels > 1 && adef->useMultipleModel);

      adaptRdataToSecondary();
    }

  if(adef->useMultipleModel)    
    category      = tr->model;
  
  index    = cdta->alias;
  data     = rdta->y;

  index[0] = -1;

#ifdef DEBUG_MSG	
  // before 
  std::cout << "aln before" << std::endl; 
  for(nat i = 1; i < _numTax + 1 ; ++i)
    {
      for(nat j = 1; j < _numSites + 1  ;  ++j)
	{
	  std::cout << inverseMeaningDNA[rdta->y[i][j]] ; 
	}
      std::cout << std::endl; 
    }

  for(nat i = 0; i < _numSites + 1 ; ++i)
    std::cout << index[i] << ","; 
  std::cout << "\n"; 
#endif

#ifdef OLD_SORT 
  // TODO i am pretty sure that this is the part that drags down the
  // parsing speed
  if(adef->compressPatterns)
    {
      for (gap = n / 2; gap > 0; gap /= 2)
	{
	  for (i = gap + 1; i <= n; i++)
	    {
	      j = i - gap;

	      do
		{
		  jj = index[j];
		  jg = index[j+gap];
		  if(adef->useMultipleModel)
		    {		     		      
		      assert(category[jj] != -1 &&
		  	     category[jg] != -1);
		     
		      flip = (category[jj] > category[jg]);
		      tied = (category[jj] == category[jg]);		     

		    }
		  else
		    {
		      flip = 0;
		      tied = 1;
		    }

		  for (k = 1; (k <= nsp) && tied; k++)
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }

		  if (flip)
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		}
	      while (flip && (j > 0));
	    }
	}
    }
#else 
    if(adef->compressPatterns)
      {
	std::sort(index  + 1 , index + _numSites + 1 ,  
		  [&](const  int &a , const int& b) -> bool 
		  {
		    // if( not ( a <= int(_numSites + 1) && b <= int(_numSites + 1) )) 
		    //   std::cout << a << "," << b << std::endl; 

		    if(adef->useMultipleModel)
		      {
			if(category[a] != category[b])
			  return category[a] < category[b]; 
		      }
		    
		    for(int i = 1; i < int(_numTax + 1) ; ++i)
		      {
			if( data[i][a] != data[i][b] ) 
			  return data[i][a] < data[i][b]; 
		      }
		    return false; 
		  });
      }
#endif


#ifdef DEBUG_MSG 
  std::cout << "aln before" << std::endl; 
  for(nat i = 1; i < _numTax + 1 ; ++i)
    {
      for(nat j = 1; j < _numSites + 1  ;  ++j)
	{
	  std::cout << inverseMeaningDNA[rdta->y[i][index[j]]] ; 
	}
      std::cout << std::endl; 
    }

  for(nat i = 0; i < _numSites + 1 ; ++i)
    std::cout << index[i] << ","; 
  std::cout << "\n"; 
#endif
}


void PhylipParser::sitecombcrunch ()
{
  int  i, sitei, j, sitej, k;
  boolean  tied;
  int 
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL;

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * (_numSites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * (_numSites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;

  if(adef->mode == PER_SITE_LL)
    {
      int i;

      assert(0);

      tr->patternPosition = (int*)malloc(sizeof(int) * _numSites);
      tr->columnPosition  = (int*)malloc(sizeof(int) * _numSites);

      for( i = 0; i < int(_numSites); i++)
	{
	  tr->patternPosition[i] = -1;
	  tr->columnPosition[i]  = -1;
	}
    }

  i = 0;
  for (j = 1; j <= int(_numSites); j++)
    {
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];
      if(!adef->compressPatterns)
	tied = 0;
      else
	{
	  if(adef->useMultipleModel)
	    {	     
	      tied = (tr->model[sitei] == tr->model[sitej]);
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);
	    }
	  else
	    tied = 1;
	}

      for (k = 1; tied && (k <= int(_numTax)); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);

      if (tied)
	{
	  if(adef->mode == PER_SITE_LL)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d from column %d also at site %d\n", i, sitei, sitej);*/
	    }


	  cdta->aliaswgt[i] += rdta->wgt[sitej];
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if (cdta->aliaswgt[i] > 0) i++;

	  if(adef->mode == PER_SITE_LL)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d is from cloumn %d\n", i, sitej);*/
	    }

	  cdta->aliaswgt[i] = rdta->wgt[sitej];
	  cdta->alias[i] = sitej;
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
    }

  cdta->endsite = i;
  if (cdta->aliaswgt[i] > 0) cdta->endsite++;

  if(adef->mode == PER_SITE_LL)
    {
      assert(0);

      for(i = 0; i < int(_numSites); i++)
	{
	  int p  = tr->patternPosition[i];
	  int c  = tr->columnPosition[i];

	  assert(p >= 0 && p < cdta->endsite);
	  assert(c >= 1 && c <= int(_numSites));
	}
    }


  if(adef->useMultipleModel)
    {
      for(i = 0; i <= int(_numSites); i++)
	{
	  tr->model[i]      = aliasModel[i];
	  tr->dataVector[i] = aliasSuperModel[i];
	}
    }

  if(adef->useMultipleModel)
    {
      free(aliasModel);
      free(aliasSuperModel);
    }     
}

boolean PhylipParser::makeweights ()
{
  int  i;

  for (i = 1; i <= int(_numSites); i++)
    cdta->alias[i] = i;

  sitesort();
  sitecombcrunch();  

  return TRUE;
}


boolean PhylipParser::makevalues()
{
  int  
    i, 
    j, 
    model, 
    modelCounter;

  unsigned char
    *y    = (unsigned char *)malloc(((size_t)_numTax) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  {
    for (i = 1; i <= int(_numTax); i++)
      for (j = 0; j < cdta->endsite; j++)   
	y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];

    free(rdta->y0);
    free(rdta->y);
   
  }

  rdta->y0 = y;
 
  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {
      tr->partitionData[0].lower = 0;

      model        = tr->model[0];
      modelCounter = 0;
     
      i            = 1;

      while(i <  cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {
	      tr->partitionData[modelCounter].upper     = (size_t)i;
	      tr->partitionData[modelCounter + 1].lower = (size_t)i;

	      model = tr->model[i];	     
	      modelCounter++;
	    }
	  i++;
	}


      tr->partitionData[tr->NumberOfModels - 1].upper = (size_t)cdta->endsite;      
    
      for(i = 0; i < tr->NumberOfModels; i++)		  
	tr->partitionData[i].width      = tr->partitionData[i].upper -  tr->partitionData[i].lower;
	 
      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;
      i            = 1;
	
      while(i < cdta->endsite)
	{	 
	  if(tr->model[i] != model)
	    {
	      model = tr->model[i];
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      
    }
  else
    {
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = (size_t)cdta->endsite;
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }

  tr->rdta       = rdta;
  tr->cdta       = cdta; 

  tr->originalCrunchedLength = tr->cdta->endsite;
    
  for(i = 0; i < int(_numTax); i++)
    tr->yVector[i + 1] = &(rdta->y0[(tr->originalCrunchedLength) * ((size_t)i)]);

  return TRUE;
}




hashNumberType  PhylipParser::hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}


void PhylipParser::addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}


stringHashtable* PhylipParser::initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

     static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
     32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
     4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
     268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}


void PhylipParser::getinput()
{
  int i;

  FILE *INFILE = fopen(alnFile.c_str(), "rb");
  
  parseHeader(INFILE); 

  tr->mxtips = _numTax;
  
  rdta->wgt             = (int *)    malloc((_numSites + 1) * sizeof(int));
  cdta->alias           = (int *)    malloc((_numSites + 1) * sizeof(int));
  cdta->aliaswgt        = (int *)    malloc((_numSites + 1) * sizeof(int)); 
  tr->model             = (int *)    calloc((_numSites + 1), sizeof(int));
  tr->initialDataVector  = (int *)    malloc((_numSites + 1) * sizeof(int));
  tr->extendedDataVector = (int *)    malloc((_numSites + 1) * sizeof(int));         
  

  assert(!adef->useWeightFile); 

  for (i = 1; i <= int(_numSites); i++)
    rdta->wgt[i] = 1;


  assert(adef->useMultipleModel); 
  
  int ref;
  
  if(_haveModelFile)
    parsePartitions();
  else 
    parseSinglePartition(modelFile);

  for(i = 1; i <= int(_numSites); i++)
    {
      ref = tr->model[i];
      tr->initialDataVector[i] = tr->initialPartitionData[ref].dataType;
    }

  assert(not adef->useSecondaryStructure); 

  tr->dataVector    = tr->initialDataVector;
  tr->partitionData = tr->initialPartitionData;

  getyspace();

  setupTree();

  if(!getdata(INFILE))
    {
      printf("Problem reading alignment file \n");
      myExit(1);
    }
      
  tr->nameHash = initStringHashTable(10 * tr->mxtips);
  for(i = 1; i <= tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  fclose(INFILE);
}


void PhylipParser::compressDNA(int *informative)
{
  size_t
    totalNodes,
    model;
   
  totalNodes = 2 * (size_t)tr->mxtips;

  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = (size_t)tr->partitionData[model].states,
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)malloc(states * sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)malloc(states * sizeof(parsimonyNumber));
      
      for(nat i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->cdta->aliaswgt[i];     

      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      size_t numByte = compressedEntriesPadded * states * totalNodes; 
      tr->partitionData[model].parsVect = (parsimonyNumber *)malloc(numByte * sizeof(parsimonyNumber));

      for(nat i = 0; i < numByte; i++)      
	tr->partitionData[model].parsVect[i] = 0;

      for(nat i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(tr->partitionData[model].parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const unsigned int 
		    *bitValue = getBitVector(tr->partitionData[model].dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->cdta->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & myMask32[k])
			    compressedValues[k] |= myMask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= myMask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	} 

      tr->partitionData[model].parsimonyLength = compressedEntriesPadded;

      free(compressedTips);
      free(compressedValues);
    }

  // formerly aligned! 
  // tr->parsimonyScore = (unsigned int*)malloc_aligned(sizeof(unsigned int) * totalNodes * tr->NumberOfModels);  
  // tr->parsimonyScore = (unsigned int*)malloc(sizeof(unsigned int) * totalNodes * tr->NumberOfModels);  
          
  // for(nat i = 0; i < totalNodes * tr->NumberOfModels; i++) 
  //   tr->parsimonyScore[i] = 0;
}

boolean PhylipParser::isInformative(int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}

void PhylipParser::determineUninformativeSites( int *informative)
{
  int 
    model,
    number = 0; 


  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */


  for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(int i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	  if(isInformative( tr->partitionData[model].dataType, i))
	    informative[i] = 1;
	  else
	    {
	      informative[i] = 0;
	      number++;
	    }  
	}      
    }

  /* printf("Uninformative Patterns: %d\n", number); */
}


Bipartition PhylipParser::allocateParsimonyDataStructures()
{
  auto informative = std::vector<int>(tr->originalCrunchedLength, 0);
 
  determineUninformativeSites( informative.data());

  compressDNA( informative.data());
 
  tr->ti = (int*)malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  
  
  auto result = Bipartition{} ;
  result.reserve(tr->originalCrunchedLength);
  for(decltype(tr->originalCrunchedLength) i = 0; i < tr->originalCrunchedLength; ++i)
    {
      if(informative[i] > 0 )
	result.set(i); 
    }

  return result;
}


void PhylipParser::writeWeights(std::ofstream &out)
{
  auto elem =  *(std::max_element(tr->cdta->aliaswgt, tr->cdta->aliaswgt + tr->originalCrunchedLength)); 

  // sorry, hard coding here 

  if(elem < std::numeric_limits<uint8_t>::max())
    {
      int len = sizeof(uint8_t); 
      myWrite(out, &len,1); 
      for(auto i = 0ull; i < tr->originalCrunchedLength; ++i)
	{
	  uint8_t val = tr->cdta->aliaswgt[i]; 
	  myWrite(out, &val, 1); 
	}
    }
  else if(elem < std::numeric_limits<uint16_t>::max())
    {
      int len = sizeof(uint16_t); 
      myWrite(out, &len,1); 
      for(auto i = 0ull; i < tr->originalCrunchedLength; ++i)
	{
	  uint16_t val = tr->cdta->aliaswgt[i]; 
	  myWrite(out, &val, 1); 
	}
    }
  else if(elem < std::numeric_limits<int32_t>::max())
    {
      int len = sizeof(uint32_t); 
      myWrite(out, &len,1); 
      for(auto i = 0ull; i < tr->originalCrunchedLength; ++i)
	{
	  uint32_t val = tr->cdta->aliaswgt[i]; 
	  myWrite(out, &val, 1); 
	}
    }
  else 
    {
      assert(0); 
      // cannot do that 
      int len = sizeof(int64_t); 
      myWrite(out, &len,1); 
      for(auto  i = 0ull; i < tr->originalCrunchedLength; ++i)
	{
	  int64_t val = tr->cdta->aliaswgt[i]; 
	  myWrite(out, &val, 1); 
	}
    }
}




void PhylipParser::writeToFile(std::string fileName) 
{
  size_t 
    i,
    model;

  auto &&out = std::ofstream(fileName, std::ios::binary); 
  
  auto fileId = std::string{ "BINARY"} ; 
  myWrite(out, fileId.c_str(), fileId.size() ); 

  myWrite(out, &(tr->mxtips), 1); 
  myWrite(out, &(tr->NumberOfModels), 1); 
  myWrite(out,&tr->originalCrunchedLength, 1) ; 

  writeWeights(out); 

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int len = strlen(tr->nameList[i]) + 1;
      myWrite(out, &len, 1); 
      myWrite(out, tr->nameList[i], len); 
    } 

  myWrite(out, tr->partitionContributions, tr->NumberOfModels); 

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      int 
	len;
	
      pInfo 
	*p = &(tr->partitionData[model]);

      myWrite(out, &(p->states),1); 
      myWrite(out, &(p->maxTipStates),1); 
      myWrite(out, &(p->lower),1); 
      myWrite(out, &(p->upper),1); 
      myWrite(out, &(p->width),1); 
      myWrite(out, &(p->dataType),1); 
      myWrite(out, &(p->protModels),1); 
      myWrite(out, &(p->protFreqs),1); 
      myWrite(out, &(p->nonGTR),1); 

      len = strlen(p->partitionName) + 1;

      myWrite(out, &len,1); 
      myWrite(out, p->partitionName, len);
    } 

#ifdef OLD_ALN_LAYOUT
  auto iter = rdta->y0; 
  for(int i = 0 ;i < tr->mxtips; ++i)
    {
      // necessary because of an overflow with a dataset consisting of
      // 200 taxa and 1e8 bp
      myWrite(out, iter  , tr->originalCrunchedLength); 
      iter += tr->originalCrunchedLength; 
    }
#else 
  auto iter = rdta->y0; 
  auto numPat = tr->originalCrunchedLength; 
  auto pattern = std::vector<uint8_t>(tr->mxtips,0); 
  for(uint64_t i = 0; i < numPat; ++i)
    {
      for(int j = 0; j < tr->mxtips; ++j)
	pattern[j] = *(iter + numPat * j + i); 
      myWrite(out, pattern.data() ,  tr->mxtips); 
    }
  
#endif

  // also write infoness 
  auto handle = _infoness.getRawBip(); 
  
  // std::cout << "INFO: again: "; 
  // for(nat i = 0 ;i < tr->originalCrunchedLength; ++i)
  //   {
  //     if(_infoness.isSet(i))
  // 	std::cout << 1; 
  //     else 
  // 	std::cout << 0 ; 
  //   }
  // std::cout << std::endl; 
  
  myWrite(out, handle.data(), handle.size()); 
  
  // for(model = 0; model < (size_t) tr->NumberOfModels; ++model)
  //   {
  //     myWrite(out, &(tr->partitionData[model].parsimonyLength),1); 
  //     size_t numBytes =tr->partitionData[model].parsimonyLength * tr->partitionData[model].states * 2 * tr->mxtips ; 
  //     myWrite(out, tr->partitionData[model].parsVect, numBytes); 
  //   }
  out.close(); 
}


void PhylipParser::parse()
{
  auto nowTime =  getTimePoint(); 

  int model;

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
  tr   = (tree *)malloc(sizeof(tree));

  defaultInit();

  getinput();  
  // auto dur = getDuration(nowTime);
  nowTime = getTimePoint(); 
  // std::cout << SOME_FIXED_PRECISION << "[" << dur << " s] got input" << std::endl; 
  makeweights();         
  // dur = getDuration(nowTime); 
  nowTime = getTimePoint(); 
  // std::cout << "[" << dur << " s] made weights" << std::endl; 
  makevalues();         
  // dur = getDuration(nowTime); 
  nowTime = getTimePoint(); 
  // std::cout << "[" << dur << " s] made values" << std::endl; 

  for(model = 0; model < tr->NumberOfModels; model++)
    {	
      int 
	states = -1,
	maxTipStates = getUndetermined(tr->partitionData[model].dataType) + 1;  	      

      switch(tr->partitionData[model].dataType)
	{
	case DNA_DATA:
	  states = getStates(tr->partitionData[model].dataType);	 
	case AA_DATA:	
	  states = getStates(tr->partitionData[model].dataType);	 
	  break;	
	default:
	  assert(0);
	}

      tr->partitionData[model].states       = states;
      tr->partitionData[model].maxTipStates = maxTipStates;
    }   


  // create the partition contributions 
  tr->partitionContributions = (double*)calloc(tr->NumberOfModels, sizeof(double)); 
  double total = 0; 
  for(int i = 0; i < tr->NumberOfModels; ++i)
    {
      double contribution = 0 ; 
      auto &partition = tr->partitionData[i]; 
      for(int j = partition.lower; j < partition.upper ;++j)
	contribution += tr->cdta->aliaswgt[j] ; 
      tr->partitionContributions[i] = contribution; 
      total += contribution; 
    }
  for(int i = 0; i < tr->NumberOfModels; ++i)
    tr->partitionContributions[i] /= total; 

  _infoness = allocateParsimonyDataStructures()  ; 

  // for(nat i = 0; i < tr->originalCrunchedLength; ++i)
  //   {
  //   if(_infoness.isSet(i))
  //     std::cout  << "1"; 
  //   else 
  //     std::cout << "0" ; 
  //   }
  // std::cout << std::endl; 
  
}
