#include "TreeInitializer.hpp"
#include <cstring>
#include <cassert>
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"

extern "C"
{
  void initializePartitions(tree *tr, FILE *byteFile); 
  void multiprocessorScheduling(tree *tr, int tid); 
  void computeFraction(tree *tr, int tid, int n); 
  void computeFractionMany(tree *tr, int tid); 
}


template<typename T>
static void byteRead(std::ifstream& in, T* result, nat num)
{
  in.read((char*)result, sizeof(T) * num ); 
}


void TreeInitializer::parseMagicNumber(std::ifstream& in)
{
  // parse the "magic number"
  char tmp[7]; 
  memset(&tmp, 0 , sizeof(char) * 7 ); 
  in.read(tmp, 6 * sizeof(char)); 
  if( std::string(tmp).compare("BINARY") != 0 )
    {
      tout << "expected different file start (maybe reparse the binary file). Got >"  << tmp << "<"<< std::endl; 
      assert(0); 
    }
}


TreeInitializer::TreeInitializer()
{
#if HAVE_PLL == 0
  adef.max_rearrange          = 21;
  adef.stepwidth              = 5;
  adef.initial                = 10;
  adef.bestTrav               = 10;
  adef.initialSet             = FALSE; 
  adef.mode                   = BIG_RAPID_MODE; 
  adef.likelihoodEpsilon      = 0.1; 
  adef.permuteTreeoptimize    = FALSE; 
  adef.perGeneBranchLengths   = TRUE;   
  adef.useCheckpoint          = FALSE;
#endif

  mask32[0] = 1 ; 
  for(nat i = 1; i < 32; ++i)
    mask32[i] = mask32[i-1] << 1 ; 
}



void TreeInitializer::unifiedInitializePartitions(TreeAln &traln, std::string byteFileName)
{
#if HAVE_PLL != 0
  initializeTreePLL(traln,byteFileName);
#else 
  initializeTreeExaML(traln, byteFileName);
#endif  
}



void TreeInitializer::initializeParsimonyVectors(TreeAln &traln, std::ifstream& byteStream)
{
  auto tr = traln.getTr(); 
  nat totalNodes = 2 * tr->mxtips; 
  nat numPart = traln.getNumberOfPartitions(); 

  for(nat model = 0; model < numPart; model++)
    {
      auto partition = traln.getPartition(model); 

      byteRead(byteStream, &(partition->parsimonyLength), 1); 
      int numBytes = totalNodes * partition->states * partition->parsimonyLength; 
      partition->parsVect = (parsimonyNumber*)exa_malloc_aligned( numBytes * sizeof(parsimonyNumber));
      memset(partition->parsVect, 0 , sizeof(parsimonyNumber) * numBytes) ;
      byteRead(byteStream, partition->parsVect, numBytes); 
    }

  for(int i = 1 ; i < 2 * tr->mxtips; ++i)
    {
      nodeptr p = tr->nodep[i]; 
      p->xPars = 1 ; 
      p->next->xPars = 0; 
      p->next->next->xPars = 0; 
    }
  

  tr->parsimonyScore = (unsigned int*)exa_malloc_aligned(sizeof(unsigned int) * totalNodes * numPart);  
  memset(tr->parsimonyScore, 0, totalNodes * numPart * sizeof(unsigned int )); 
  
}


void TreeInitializer::readPartitions(TreeAln & traln, std::ifstream &byteStream)
{
  
  for(nat model = 0; model < traln.getNumberOfPartitions() ; model++)
    {      
      int 
	len = 0;

      auto p = traln.getPartition(model); 
      
      byteRead(byteStream, &p->states, 1); 
      byteRead(byteStream, &p->maxTipStates, 1); 
      byteRead(byteStream, &p->lower, 1); 
      byteRead(byteStream, &p->upper, 1); 
      byteRead(byteStream, &p->width, 1); 
      byteRead(byteStream, &p->dataType, 1); 
      byteRead(byteStream, &p->protModels, 1); 
      byteRead(byteStream, &p->protFreqs, 1); 
      byteRead(byteStream, &p->nonGTR, 1); 
      
      byteRead(byteStream, &len, 1); 
      p->partitionName = (char*)malloc(sizeof(char) * len);
      byteRead(byteStream, p->partitionName, len); 
    }
}



#if HAVE_PLL == 0
void TreeInitializer::initializePartitionsExaml(TreeAln &traln, std::ifstream &byteStream)
{ 
  auto tr = traln.getTr(); 

  size_t
    i,
    j,
    width,
    model,
    countOffset,
    myLength = 0;

  int
    maxCategories;

  unsigned char 
    *y;

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    tr->partitionData[model].width = 0;

  if(tr->manyPartitions)
    {
      multiprocessorScheduling(tr, processID);  
      computeFractionMany(tr, processID); 
    }
  else
    computeFraction(tr, processID, processes);
  	   
  maxCategories = tr->maxCategories;

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {                       
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[model])); 

      width = tr->partitionData[model].width;
	
      /* 
	 globalScaler needs to be 2 * tr->mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
       */

      tr->partitionData[model].globalScaler    = (unsigned int *)calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      tr->partitionData[model].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[model].EIGN              = (double*)malloc(pl->eignLength * sizeof(double));
      tr->partitionData[model].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[model].EI                = (double*)malloc(pl->eiLength * sizeof(double));
      
      tr->partitionData[model].substRates        = (double *)malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[model].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].empiricalFrequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[model].tipVector         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));

      if(tr->partitionData[model].protModels == LG4)      
	{	  	  
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {	    
	      tr->partitionData[model].EIGN_LG4[k]              = (double*)malloc(pl->eignLength * sizeof(double));
	      tr->partitionData[model].EV_LG4[k]                = (double*)malloc_aligned(pl->evLength * sizeof(double));
	      tr->partitionData[model].EI_LG4[k]                = (double*)malloc(pl->eiLength * sizeof(double));
	      tr->partitionData[model].substRates_LG4[k]        = (double *)malloc(pl->substRatesLength * sizeof(double));
	      tr->partitionData[model].frequencies_LG4[k]       = (double*)malloc(pl->frequenciesLength * sizeof(double));
	      tr->partitionData[model].tipVector_LG4[k]         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
	    }
	}


      tr->partitionData[model].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[model].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));
      
      tr->partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * tr->maxCategories);
            
      tr->partitionData[model].nonGTR = FALSE;            

      tr->partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[model].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

      
      tr->partitionData[model].xVector = (double **)malloc(sizeof(double*) * tr->mxtips);   
      	
      for(j = 0; j < (size_t)tr->mxtips; j++)	        	  	  	  	 
	  tr->partitionData[model].xVector[j]   = (double*)NULL;   

      tr->partitionData[model].xSpaceVector = (size_t *)calloc(tr->mxtips, sizeof(size_t));  

      tr->partitionData[model].sumBuffer = (double *)malloc_aligned(width *
									   (size_t)(tr->partitionData[model].states) *
									   discreteRateCategories(tr->rateHetModel) *
									   sizeof(double));
	    
      tr->partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      tr->partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));

      if(width > 0 && tr->saveMemory)
	{
	  tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
	    
	  tr->partitionData[model].gapVector = (unsigned int*)calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));	  	    	  	  
	    
	  tr->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr->mxtips) *								      
									       ((size_t)(tr->partitionData[model].states)) *
									       discreteRateCategories(tr->rateHetModel) * sizeof(double));
	}
      else
	{
	   tr->partitionData[model].gapVectorLength = 0;
	    
	   tr->partitionData[model].gapVector = (unsigned int*)NULL; 	  	    	   
	    
	   tr->partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }

        
  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    myLength += tr->partitionData[model].width;         
   
  /* assign local memory for storing sequence data */

  tr->y_ptr = (unsigned char *)malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);
   
  for(i = 0; i < (size_t)tr->mxtips; i++)
    {
      for(model = 0, countOffset = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  tr->partitionData[model].yVector[i+1]   = &tr->y_ptr[i * myLength + countOffset];
	  countOffset +=  tr->partitionData[model].width;
	}
      assert(countOffset == myLength);
    }

  /* figure in data */

  if(tr->manyPartitions)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  if(isThisMyPartition(tr, processID, model))
	    {
	      width = tr->partitionData[model].upper - tr->partitionData[model].lower;	     
	      
	      memcpy(&(tr->partitionData[model].wgt[0]), &(tr->aliaswgt[tr->partitionData[model].lower]), sizeof(int) * width);
	    }
	}
    }
  else
    {
      size_t 	   
	globalCounter, 
	r, 
	localCounter;
      
      for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
	    {
	      if(r % (size_t)processes == (size_t)processID)
		{
		  tr->partitionData[model].wgt[localCounter] = tr->aliaswgt[globalCounter]; 
		  
		  localCounter++;
		}
	      globalCounter++;
	    }
	  assert(localCounter == nat(tr->partitionData[model].width));
	}   
      assert(globalCounter == nat(tr->originalCrunchedLength));
    }
   
  /* set up the averaged frac changes per partition such that no further reading accesses to aliaswgt are necessary
     and we can free the array for the GAMMA model */

  if(tr->NumberOfModels > 1)
    {        
      size_t
	*modelWeights = (size_t *)calloc(tr->NumberOfModels, sizeof(size_t)),
	wgtsum = 0;  
      
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)      
	 {
	   size_t
	     lower = tr->partitionData[model].lower,
	     upper = tr->partitionData[model].upper,
	     i;
	   
	   for(i = lower; i < upper; i++)
	     {
	       modelWeights[model] += (size_t)tr->aliaswgt[i];
	       wgtsum              += (size_t)tr->aliaswgt[i];
	     }
	 }
       
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)      	
	 tr->partitionContributions[model] = ((double)modelWeights[model]) / ((double)wgtsum); 
       
       free(modelWeights);
    }

  if(tr->rateHetModel == GAMMA)
    {      
      free(tr->aliaswgt);
      tr->aliaswgt = NULL; 
    }


  y = (unsigned char *)malloc(sizeof(unsigned char) * tr->originalCrunchedLength);

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      byteRead(byteStream, y, tr->originalCrunchedLength );
	
      if(tr->manyPartitions)
	{
	  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	    {
	      if(isThisMyPartition(tr, processID, model))	  
		{
		  memcpy(tr->partitionData[model].yVector[i], &(y[tr->partitionData[model].lower]), sizeof(unsigned char) * tr->partitionData[model].width);					    
		  assert(tr->partitionData[model].width == tr->partitionData[model].upper - tr->partitionData[model].lower);
		}
	      else
		assert(tr->partitionData[model].width == 0);
	    }
	}
      else
	{
	  size_t 	  
	    globalCounter, 
	    r, 
	    localCounter;

	  for(model = 0, globalCounter = 0; model < (size_t)tr->NumberOfModels; model++)
	    {
	      for(localCounter = 0, r = (size_t)tr->partitionData[model].lower;  r < (size_t)tr->partitionData[model].upper; r++)
		{
		  if(r % (size_t)processes == (size_t)processID)
		    {		      
		      tr->partitionData[model].yVector[i][localCounter] = y[globalCounter]; 	     
		      
		      localCounter++;
		    }
		  globalCounter++;
		}
	      
	      assert(localCounter == nat(tr->partitionData[model].width));
	    }

	  assert(globalCounter == nat(tr->originalCrunchedLength));
	}
    }

  free(y);
    
  /* initialize gap bit vectors at tips when memory saving option is enabled */
  
  if(tr->saveMemory)
    {
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{
	  int        
	    undetermined = getUndetermined(tr->partitionData[model].dataType);
	  	 
	  width =  tr->partitionData[model].width;
	    
	  if(width > 0)
	    {	   	    	      	    	     
	      for(j = 1; j <= (size_t)(tr->mxtips); j++)
		for(i = 0; i < width; i++)
		  if(tr->partitionData[model].yVector[j][i] == undetermined)
		    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }
}

void TreeInitializer::initializeTreeExaML(TreeAln &traln, std::string byteFileName )
{
  auto tr = traln.getTr();  

  size_t 
    i ;

  auto&& byteFile = std::ifstream{byteFileName}; 

  parseMagicNumber(byteFile); 

  byteRead(byteFile, &tr->mxtips,1); // 1
  byteRead(byteFile, &tr->NumberOfModels, 1);	      // 2
  NUM_BRANCHES = tr->NumberOfModels; 
  byteRead(byteFile, &tr->gapyness,1);		      // 4
  byteRead(byteFile, &tr->originalCrunchedLength, 1); // 3

  tr->numBranches = tr->NumberOfModels;

  tr->aliaswgt                   = (int *)malloc(tr->originalCrunchedLength * sizeof(int));
  byteRead(byteFile, tr->aliaswgt, tr->originalCrunchedLength); 

  tr->executeModel   = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
  for(i = 0; i < (size_t)tr->NumberOfModels; i++)
    tr->executeModel[i] = TRUE;
   
  setupTree(tr); 		// TREE SETUP 

  for(i = 1; i <= (size_t)tr->mxtips; i++)
    {
      int 
	len = 0;

      byteRead(byteFile, &len, 1); 
      tr->nameList[i] = (char*)malloc(sizeof(char) * len);
      byteRead(byteFile, tr->nameList[i] , len); 
      addword(tr->nameList[i], tr->nameHash, i);        
    }  

  readPartitions(traln, byteFile);     

  initializePartitionsExaml(traln, byteFile);

  initializeParsimonyVectors(traln, byteFile);

  tr->ti = (int*) malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  

  unifiedModelInit(traln);
}

#else  

void TreeInitializer::initializeTreePLL(TreeAln &traln, std::string byteFileName)
{
  initializePartitionsPLL(traln, byteFileName );
  initializePartitionsSequential(traln.getTr(), traln.getPartitionsPtr());

  unifiedModelInit(traln);

  auto ptr = traln.getPartitionsPtr(); 
  ptr->perGeneBranchLengths = TRUE; 
} 


void TreeInitializer::initializePartitionsPLL(TreeAln &traln, std::string byteFileName)
{
  unsigned char *y;

  auto &&byteStream  = std::ifstream{byteFileName, std::ios::binary }; 

  parseMagicNumber(byteStream); 

  auto ptr = traln.getPartitionsPtr(); 
  auto& partitions = *ptr; 
  auto trPtr = traln.getTr(); 
  auto& tr = *trPtr; 
  
  ptr->perGeneBranchLengths = 1;
  
  byteRead(byteStream, &tr.mxtips, 1); // 1
  byteRead(byteStream, &partitions.numberOfPartitions, 1) ; // 2
  NUM_BRANCHES = partitions.numberOfPartitions; 
  byteRead(byteStream, &tr.gapyness, 1); // 3
  byteRead(byteStream, &tr.originalCrunchedLength,1); // 4

  tr.aliaswgt = (int *)exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));
  byteRead(byteStream, tr.aliaswgt, tr.originalCrunchedLength); // 5

  partitions.partitionData = (pInfo**)exa_calloc(partitions.numberOfPartitions,  sizeof(pInfo*));

  tr.rateCategory    = (int *)    exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));	  
  tr.patrat          = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double));
  tr.patratStored    = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double)); 
  tr.lhs             = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double)); 

  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr.originalCrunchedLength) * ((size_t)tr.mxtips));
  tr.yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr.mxtips + 1)));

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    tr.yVector[i] = &y[(i - 1) *  (size_t)tr.originalCrunchedLength]; 

  setupTree(&tr, FALSE, &partitions); // TREE SETUP 

  for(int i = 0; i < partitions.numberOfPartitions; i++)
    partitions.partitionData[i]->executeModel = TRUE;

  /* data structures for convergence criterion need to be initialized after! setupTree */
  
  tr.nameList = (char**)exa_calloc(2 * tr.mxtips , sizeof(char*)); 

  for( int i = 1; i <= tr.mxtips; i++)
    {
      int len = 0;
      byteRead(byteStream, &len, 1); 

      tr.nameList[i] = (char*)exa_calloc(len, sizeof(char));
      byteRead(byteStream, tr.nameList[i], len); 
    }  

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    addword(tr.nameList[i], tr.nameHash, i);
 

  readPartitions(traln, byteStream);

  byteRead(byteStream, y, tr.originalCrunchedLength * tr.mxtips); 
}

#endif


void TreeInitializer::unifiedModelInit(TreeAln &traln)
{
  double** empFreqs = new double*[traln.getNumberOfPartitions()]; 
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
      auto partition = traln.getPartition(i); 
      empFreqs[i] = new double[partition->states]; 
      for(int j = 0; j < partition->states ; ++j)
	empFreqs[i][j] = 1. / partition->states; 
    }

#if HAVE_PLL == 0 
  initModel(traln.getTr(), empFreqs);
#else 
  initModel(traln.getTr(), empFreqs, traln.getPartitionsPtr());
#endif
  
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    delete [] empFreqs[i];
}
