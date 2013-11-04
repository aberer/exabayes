#include "TreeInitializer.hpp"

#include <cstring>
#include <cassert>
#include "TreeAln.hpp"
#include "GlobalVariables.hpp"

extern "C"
{
  void multiprocessorScheduling(tree *tr, int tid); 
  void computeFraction(tree *tr, int tid, int n); 
  void computeFractionMany(tree *tr, int tid); 
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


#if HAVE_PLL == 0
void TreeInitializer::distributeData(TreeAln &traln)
{
  auto &tr = traln.getTrHandle();

  /* figure in data */
  if(tr.manyPartitions)
    {
      for(int model = 0; model < tr.NumberOfModels; model++)
	{
	  if(isThisMyPartition(&tr, processID, model))
	    {
	      int width = tr.partitionData[model].upper - tr.partitionData[model].lower;	     

	      memcpy(&(tr.partitionData[model].wgt[0]), &(tr.aliaswgt[tr.partitionData[model].lower]), sizeof(int) * width);
	    }
	}
    }
  else
    {
      size_t 	   
	globalCounter = 0, 
	r, 
	localCounter;
      
       for(int model = 0; model < tr.NumberOfModels; model++)
	{
	  for(localCounter = 0, r = (size_t)tr.partitionData[model].lower;  r < (size_t)tr.partitionData[model].upper; r++)
	    {
	      if(r % (size_t)processes == (size_t)processID)
		{
		  tr.partitionData[model].wgt[localCounter] = tr.aliaswgt[globalCounter]; 
		  
		  localCounter++;
		}
	      globalCounter++;
	    }
	  assert(localCounter == nat(tr.partitionData[model].width));
	}   
      assert(globalCounter == nat(tr.originalCrunchedLength));
    }

  unsigned char *y = (unsigned char *)malloc(sizeof(unsigned char) * tr.originalCrunchedLength);

  nat ctr = 0; 
  for(int i = 1; i <= tr.mxtips; i++)
    {
      _initResPtr->fillAlnPart(y,tr.originalCrunchedLength,ctr);
	
      if(tr.manyPartitions)
	{
	  for(int model = 0; model < tr.NumberOfModels; model++)
	    {
	      if(isThisMyPartition(&tr, processID, model))	  
		{
		  memcpy(tr.partitionData[model].yVector[i], &(y[tr.partitionData[model].lower]), sizeof(unsigned char) * tr.partitionData[model].width);					    
		  assert(tr.partitionData[model].width == tr.partitionData[model].upper - tr.partitionData[model].lower);
		}
	      else
		assert(tr.partitionData[model].width == 0);
	    }
	}
      else
	{
	  size_t 	  
	    globalCounter = 0, 
	    r, 
	    localCounter;

	  for(int model = 0 ; model < tr.NumberOfModels; model++)
	    {
	      for(localCounter = 0, r = (size_t)tr.partitionData[model].lower;  r < (size_t)tr.partitionData[model].upper; r++)
		{
		  if(r % (size_t)processes == (size_t)processID)
		    {		      
		      tr.partitionData[model].yVector[i][localCounter] = y[globalCounter]; 	     
		      
		      localCounter++;
		    }
		  globalCounter++;
		}
	      
	      assert(localCounter == nat(tr.partitionData[model].width));
	    }

	  assert(globalCounter == nat(tr.originalCrunchedLength));
	}
    }
  free(y);
}


#else 

void TreeInitializer::distributeData(TreeAln &traln)
{
  auto &tr = traln.getTrHandle(); 

  for(nat model = 0; model < traln.getNumberOfPartitions(); model++)
    {
      auto &partition = traln.getPartition(model);
      size_t
  	j;
      size_t lower = partition.lower;
      size_t width = partition.upper - lower;

      for(j = 1; j <= (size_t)tr.mxtips; j++)
	{
	  // partition.yVector[j] = (unsigned char*)exa_calloc(width , sizeof(unsigned char)); 
	  auto start = tr.yVector[j] + lower   ; 
	  std::copy(start , start + width , partition.yVector[j]);
	}

      memcpy((void*)(partition.wgt),         (void*)(&(tr.aliaswgt[lower])),      sizeof(int) * width);
    }
}

#endif


void TreeInitializer::copyWeightsAndAln(TreeAln &traln )
{
  // assert(0); 
  // TODO initialize memory first  

  _initResPtr->initWeightsAndAln(traln);
}



void TreeInitializer::determinePartitionContribution(TreeAln &traln)
{
  auto &tr = traln.getTrHandle();

  auto localSums = std::vector<nat>{}; 
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i )
    {
      auto& partition = traln.getPartition(i); 
      auto localSum = nat{0} ;
      for(nat j = partition.lower; j < nat(partition.upper); ++j)
	localSum += tr.aliaswgt[j]; 
      localSums.push_back(localSum); 
    }    
  auto totalSum = std::accumulate(localSums.begin(), localSums.end(), 0.);
  
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
#if HAVE_PLL == 0
      tr.partitionContributions[i] = double(localSums[i])  / double(totalSum); 
#else 
      auto &partition = traln.getPartition(i); 
      partition.partitionContribution = double(localSums[i])  / double(totalSum);
#endif
    }
}






TreeInitializer::TreeInitializer(std::unique_ptr<InitializationResource> initRes)
  : _initResPtr(std::move(initRes))
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


void TreeInitializer::unifiedInitializePartitions(TreeAln &traln)
{
#if HAVE_PLL != 0
  initializeTreePLL(traln);
#else 
  initializeTreeExaML(traln);
#endif  
}


void TreeInitializer::initializeParsimonyVectors(TreeAln &traln)
{ 
  auto& tr = traln.getTrHandle(); 
  nat totalNodes = 2 * tr.mxtips; 
  nat numPart = traln.getNumberOfPartitions(); 

  nat ctr = 0; 
  for(nat model = 0; model < numPart; model++)
    {
      auto& partition = traln.getPartition(model); 
      _initResPtr->fillParsVect(partition.parsVect, partition.parsimonyLength, totalNodes * partition.states, ctr);
    }

  for(int i = 1 ; i < 2 * tr.mxtips; ++i)
    {
      nodeptr p = tr.nodep[i]; 
      p->xPars = 1 ; 
      p->next->xPars = 0; 
      p->next->next->xPars = 0; 
    }

  tr.parsimonyScore = (unsigned int*)exa_malloc_aligned(sizeof(unsigned int) * totalNodes * numPart);  
  memset(tr.parsimonyScore, 0, totalNodes * numPart * sizeof(unsigned int )); 
  tr.ti = (int*) malloc(sizeof(int) * 4 * (size_t)tr.mxtips);  
}


#if HAVE_PLL == 0
void TreeInitializer::initializePartitionsExaml(TreeAln &traln)
{ 
  auto& tr = traln.getTrHandle(); 

  for(int model = 0; model < tr.NumberOfModels; model++)
    tr.partitionData[model].width = 0;

  if(tr.manyPartitions)
    {
      multiprocessorScheduling(&tr, processID);  
      computeFractionMany(&tr, processID); 
    }
  else
    computeFraction(&tr, processID, processes);

  int maxCategories = tr.maxCategories;

  for(int model = 0; model < tr.NumberOfModels; model++)
    {                       
      const partitionLengths 
	*pl = getPartitionLengths(&(tr.partitionData[model])); 

      int width = tr.partitionData[model].width;

      /* 
	 globalScaler needs to be 2 * tr.mxtips such that scalers of inner AND tip nodes can be added without a case switch
	 to this end, it must also be initialized with zeros -> calloc
      */

      tr.partitionData[model].globalScaler    = (unsigned int *)calloc(2 * tr.mxtips, sizeof(unsigned int));  	         

      tr.partitionData[model].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr.partitionData[model].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr.partitionData[model].EIGN              = (double*)exa_calloc(pl->eignLength , sizeof(double));
      tr.partitionData[model].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr.partitionData[model].EI                = (double*)exa_calloc(pl->eiLength,  sizeof(double));

      tr.partitionData[model].substRates        = (double *)exa_calloc(pl->substRatesLength, sizeof(double));
      tr.partitionData[model].frequencies       = (double*)exa_calloc(pl->frequenciesLength, sizeof(double));
      tr.partitionData[model].empiricalFrequencies       = (double*)exa_calloc(pl->frequenciesLength, sizeof(double));
      tr.partitionData[model].tipVector         = (double *)exa_malloc_aligned(pl->tipVectorLength * sizeof(double));

      if(tr.partitionData[model].protModels == LG4)      
	{	  	  
	  int 
	    k;

	  for(k = 0; k < 4; k++)
	    {	    
	      tr.partitionData[model].EIGN_LG4[k]              = (double*)exa_calloc(pl->eignLength,  sizeof(double));
	      tr.partitionData[model].EV_LG4[k]                = (double*)malloc_aligned(pl->evLength * sizeof(double));
	      tr.partitionData[model].EI_LG4[k]                = (double*)exa_calloc(pl->eiLength, sizeof(double));
	      tr.partitionData[model].substRates_LG4[k]        = (double *)exa_calloc(pl->substRatesLength, sizeof(double));
	      tr.partitionData[model].frequencies_LG4[k]       = (double*)exa_calloc(pl->frequenciesLength, sizeof(double));
	      tr.partitionData[model].tipVector_LG4[k]         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
	    }
	}

      tr.partitionData[model].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr.partitionData[model].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));

      tr.partitionData[model].perSiteRates      = (double *)malloc(sizeof(double) * tr.maxCategories);

      tr.partitionData[model].nonGTR = FALSE;            

      tr.partitionData[model].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr.partitionData[model].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr.mxtips + 1));


      tr.partitionData[model].xVector = (double **)malloc(sizeof(double*) * tr.mxtips);   

      for(int j = 0; j < tr.mxtips; j++)	        	  	  	  	 
	tr.partitionData[model].xVector[j]   = (double*)NULL;   

      tr.partitionData[model].xSpaceVector = (size_t *)calloc(tr.mxtips, sizeof(size_t));  

      tr.partitionData[model].sumBuffer = (double *)malloc_aligned(width *
								   (size_t)(tr.partitionData[model].states) *
								   discreteRateCategories(tr.rateHetModel) *
								   sizeof(double));

      tr.partitionData[model].wgt = (int *)malloc_aligned(width * sizeof(int));	  

      /* rateCategory must be assigned using calloc() at start up there is only one rate category 0 for all sites */

      tr.partitionData[model].rateCategory = (int *)calloc(width, sizeof(int));

      if(width > 0 && tr.saveMemory)
	{
	  tr.partitionData[model].gapVectorLength = ((int)width / 32) + 1;

	  tr.partitionData[model].gapVector = (unsigned int*)calloc(tr.partitionData[model].gapVectorLength * 2 * tr.mxtips, sizeof(unsigned int));	  	    	  	  

	  tr.partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr.mxtips) *								      
								       ((size_t)(tr.partitionData[model].states)) *
								       discreteRateCategories(tr.rateHetModel) * sizeof(double));
	}
      else
	{
	  tr.partitionData[model].gapVectorLength = 0;

	  tr.partitionData[model].gapVector = (unsigned int*)NULL; 	  	    	   

	  tr.partitionData[model].gapColumn = (double*)NULL;	    	    	   
	}              
    }


  nat myLength = 0 ;
  for(int model = 0; model < tr.NumberOfModels; model++)
    myLength += tr.partitionData[model].width;         

  /* assign local memory for storing sequence data */

  // reserve the memory for the alignment   
  for(int i = 0; i < tr.mxtips; i++)
    {
      for(int model = 0 ; model < tr.NumberOfModels; model++)
	tr.partitionData[model].yVector[i+1] = (unsigned char*) exa_calloc(tr.partitionData[model].width, sizeof(unsigned char));
    }

  // distribute the actual alignment 
  if(_initResPtr->isDataAreDistributed())
    copyWeightsAndAln(traln); 
  else 
    distributeData(traln);

  if(tr.saveMemory)
    {
      for(int model = 0; model < tr.NumberOfModels; model++)
	{
	  int        
	    undetermined = getUndetermined(tr.partitionData[model].dataType);
	  	 
	  int width =  tr.partitionData[model].width;
	    
	  if(width > 0)
	    {	   	    	      	    	     
	      for(int j = 1; j <= tr.mxtips; j++)
		for(int i = 0; i < width; i++)
		  if(tr.partitionData[model].yVector[j][i] == undetermined)
		    tr.partitionData[model].gapVector[tr.partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];	    
	    }     
	}
    }

  if(not _initResPtr->isDataAreDistributed()) 
    {
      exa_free(tr.aliaswgt);
      tr.aliaswgt = NULL; 
    }

  if(tr.yVector)
    {
      // should not be the case 
      assert(0); 
      exa_free(tr.yVector); 
      tr.yVector = NULL; 
    }
}


void TreeInitializer::initializeTreeExaML(TreeAln &traln )
{
  auto& tr = traln.getTrHandle();  

  size_t 
    i ;

  std::tie(tr.mxtips,tr.NumberOfModels, 
	   tr.gapyness, tr.originalCrunchedLength ) = _initResPtr->getGlobalInfo();

  NUM_BRANCHES = tr.NumberOfModels; 
  tr.numBranches = tr.NumberOfModels;

  if( not _initResPtr->isDataAreDistributed())
    {
      tr.aliaswgt                   = (int *)malloc(tr.originalCrunchedLength * sizeof(int));
      _initResPtr->fillAliasWgt(tr.aliaswgt, tr.originalCrunchedLength);
    }

  tr.executeModel   = (boolean *)malloc(sizeof(boolean) * tr.NumberOfModels);
  for(i = 0; i < (size_t)tr.NumberOfModels; i++)
    tr.executeModel[i] = TRUE;

  tr.fracchanges = (double*)exa_calloc(tr.NumberOfModels, sizeof(double)); 
  for(decltype(tr.NumberOfModels) i = 0; i < tr.NumberOfModels; i++)
    tr.fracchanges[i] = -1.0;
  tr.fracchange = -1.0;

  tr.partitionData = (pInfo*)calloc(tr.NumberOfModels, sizeof(pInfo));

  tr.perPartitionLH = (double *)malloc(sizeof(double) * tr.NumberOfModels);
  for(decltype(tr.NumberOfModels) i = 0; i < tr.NumberOfModels; i++)    
    tr.perPartitionLH[i] = 0.0;	    

  traln.setTaxa(_initResPtr->getTaxonNames(tr.mxtips)); 

  auto partConts = _initResPtr->getPartitionContributions(traln.getNumberOfPartitions());
  tr.partitionContributions = (double *)malloc(sizeof(double) * tr.NumberOfModels);
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    tr.partitionContributions[i] = partConts[i]; 
  
  for(nat model = 0; model < traln.getNumberOfPartitions() ; model++)
    {      
      auto& p = traln.getPartition(model); 
      _initResPtr->fillPartition(p,model); 
    }

  initializePartitionsExaml(traln);
  initializeParsimonyVectors(traln);
  
  unifiedModelInit(traln);
}

#else  

void TreeInitializer::initializeTreePLL(TreeAln &traln)
{
  initializePartitionsPLL(traln );
  unifiedModelInit(traln);
  auto &ptr = traln.getPartitionsHandle(); 
  auto &tr =  traln.getTrHandle(); 
  ptr.perGeneBranchLengths = TRUE; 
  initializeParsimonyVectors(traln);
  
  if(tr.aliaswgt)
    {
      exa_free(tr.aliaswgt); 
      tr.aliaswgt = NULL; 
    }
  exa_free(tr.yVector); 
  tr.yVector = NULL; 
} 


void TreeInitializer::initializePartitionsPLL(TreeAln &traln)
{
  unsigned char *y;

  auto& partitions = traln.getPartitionsHandle(); 
  auto& tr = traln.getTrHandle(); 
  
  partitions.perGeneBranchLengths = 1;
  
  std::tie(tr.mxtips , partitions.numberOfPartitions, tr.gapyness, tr.originalCrunchedLength) = _initResPtr->getGlobalInfo(); 

  NUM_BRANCHES = partitions.numberOfPartitions; 

  tr.aliaswgt = (int *)exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));

  _initResPtr->fillAliasWgt(tr.aliaswgt, tr.originalCrunchedLength);

  partitions.partitionData = (pInfo**)exa_calloc(partitions.numberOfPartitions,  sizeof(pInfo*));

  tr.rateCategory    = (int *)    exa_malloc((size_t)tr.originalCrunchedLength * sizeof(int));	  
  tr.patrat          = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double));
  tr.patratStored    = (double*)  exa_malloc((size_t)tr.originalCrunchedLength * sizeof(double)); 

  y = (unsigned char *)exa_malloc(sizeof(unsigned char) * ((size_t)tr.originalCrunchedLength) * ((size_t)tr.mxtips));
  tr.yVector = (unsigned char **)exa_malloc(sizeof(unsigned char *) * ((size_t)(tr.mxtips + 1)));

  for( nat i = 1; i <= (size_t)tr.mxtips; i++)
    tr.yVector[i] = &y[(i - 1) *  (size_t)tr.originalCrunchedLength]; 

  for (nat i = 0; i < nat(partitions.numberOfPartitions); i++)
    {
      partitions.partitionData[i] = (pInfo*)exa_malloc (sizeof(pInfo));
      auto& partition = *(partitions.partitionData[i]); 
      partition.partitionContribution = -1.0;
      partition.partitionLH = 0.0;
      partition.fracchange = 1.0;
      partition.executeModel = TRUE;
    }

  traln.setTaxa(_initResPtr->getTaxonNames(tr.mxtips));

  auto partConts = _initResPtr->getPartitionContributions(traln.getNumberOfPartitions());
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
      auto& partition = traln.getPartition(i); 
      partition.partitionContribution = partConts.at(i); 
    }

  for(nat model = 0; model < traln.getNumberOfPartitions() ; model++)
    {      
      auto& p = traln.getPartition(model); 
      _initResPtr->fillPartition(p, model); 
    }

  if(not _initResPtr->isDataAreDistributed())
    {
      nat ctr = 0; 
      _initResPtr->fillAlnPart(y,tr.originalCrunchedLength * tr.mxtips, ctr); 
    }

  auto& pr = traln.getPartitionsHandle();

  for(int model = 0; model < pr.numberOfPartitions; model++)
    assert(pr.partitionData[model]->width == pr.partitionData[model]->upper - pr.partitionData[model]->lower);

  size_t 
    maxCategories = (size_t)tr.maxCategories;

  for(int model = 0; model < pr.numberOfPartitions; model++)
    {
      size_t 
	width = pr.partitionData[model]->width;

      const partitionLengths 
	*pl = getPartitionLengths(pr.partitionData[model]);
	
      auto& partition = traln.getPartition(model);

      partition.globalScaler       = (unsigned int *)exa_calloc(2 *(size_t)tr.mxtips, sizeof(unsigned int));
      partition.left              = (double *)exa_malloc_aligned((size_t)pl->leftLength * (maxCategories + 1) * sizeof(double));
      partition.right             = (double *)exa_malloc_aligned((size_t)pl->rightLength * (maxCategories + 1) * sizeof(double));
      partition.EIGN              = (double*)exa_malloc((size_t)pl->eignLength * sizeof(double));
      partition.EV                = (double*)exa_malloc_aligned((size_t)pl->evLength * sizeof(double));
      partition.EI                = (double*)exa_malloc((size_t)pl->eiLength * sizeof(double));

      partition.substRates        = (double *)exa_malloc((size_t)pl->substRatesLength * sizeof(double));
      partition.frequencies       = (double*)exa_malloc((size_t)pl->frequenciesLength * sizeof(double));
      partition.empiricalFrequencies       = (double*)exa_malloc((size_t)pl->frequenciesLength * sizeof(double));
      partition.tipVector         = (double *)exa_malloc_aligned((size_t)pl->tipVectorLength * sizeof(double));

      partition.perSiteRates      = (double *)exa_malloc(sizeof(double) * maxCategories);
      partition.nonGTR = FALSE;
      partition.gammaRates = (double*)exa_malloc(sizeof(double) * 4);
      partition.yVector = (unsigned char **)exa_malloc(sizeof(unsigned char*) * ((size_t)tr.mxtips + 1));
      partition.xVector = (double **)exa_calloc(sizeof(double*), (size_t)tr.mxtips);
      partition.xSpaceVector = (size_t *)exa_calloc((size_t)tr.mxtips, sizeof(size_t));
      partition.wgt = (int *)exa_malloc_aligned(width * sizeof(int));
      partition.rateCategory = (int *)exa_calloc(width, sizeof(int));

      if(width > 0 && tr.saveMemory)
	{
	  partition.gapVectorLength = ((int)width / 32) + 1;
	  assert(4 == sizeof(unsigned int));
	  partition.gapVector = (unsigned int*)exa_calloc((size_t)partition.gapVectorLength * 2 * (size_t)tr.mxtips, sizeof(unsigned int));
	  partition.gapColumn = (double *)exa_malloc_aligned(((size_t)tr.mxtips) * ((size_t)(partition.states)) * discreteRateCategories(tr.rateHetModel) * sizeof(double));
	}
      else
	{
	  partition.gapVectorLength = 0;
	  partition.gapVector = (unsigned int*)NULL;
	  partition.gapColumn = (double*)NULL;
	}              
    }


  for(int i = 0; i < tr.mxtips; i++)
    {
      for(nat model = 0 ; model < traln.getNumberOfPartitions(); ++model)
	{
	  auto &partition = traln.getPartition(model);
	  partition.yVector[i+1] = (unsigned char*) exa_calloc(partition.width, sizeof(unsigned char));
	}
    }

  
  if(_initResPtr->isDataAreDistributed())
    copyWeightsAndAln(traln); 
  else 
    distributeData(traln);

  
  /* initialize gap bit vectors at tips when memory saving option is enabled */
  if(tr.saveMemory)
    {
      for(int model = 0; model < pr.numberOfPartitions; model++)
	{
	  auto &partition  = traln.getPartition(model); 
	  int        
	    undetermined = getUndetermined(partition.dataType);
	  nat
	    width =  partition.width;

	  if(width > 0)
	    {	   	    	      	    	     
	      for(int j = 1; j <= tr.mxtips; j++)
		for(nat i = 0; i < width; i++)
		  if(partition.yVector[j][i] == undetermined)
		    partition.gapVector[partition.gapVectorLength * j + i / 32] |= mask32[i % 32];
	    }     
	}
    }
}

#endif


void TreeInitializer::unifiedModelInit(TreeAln &traln)
{
  double** empFreqs = new double*[traln.getNumberOfPartitions()]; 
  for(nat i = 0; i < traln.getNumberOfPartitions() ; ++i)
    {
      auto& partition = traln.getPartition(i); 
      empFreqs[i] = new double[partition.states]; 
      for(int j = 0; j < partition.states ; ++j)
	empFreqs[i][j] = 1. / partition.states; 
    }

#if HAVE_PLL == 0 
  initModel(&(traln.getTrHandle()), empFreqs);
#else 
  initModel(&(traln.getTrHandle()), empFreqs, &(traln.getPartitionsHandle()));
#endif
  
  for(nat i = 0; i < traln.getNumberOfPartitions(); ++i)
    delete [] empFreqs[i];
}



void TreeInitializer::setupTheTree(tree &tr)
{
  nodeptr  p0, p, q;

  tr.bigCutoff = FALSE;

  tr.maxCategories = MAX(4, tr.categories);

  int tips  = (size_t)tr.mxtips;
  int inter = (size_t)(tr.mxtips - 1);

  tr.treeStringLength = tr.mxtips * (nmlngth+128) + 256 + tr.mxtips * 2;

  tr.tree_string  = (char*)exa_calloc((size_t)tr.treeStringLength, sizeof(char)); 
  tr.tree0 = (char*)exa_calloc((size_t)tr.treeStringLength, sizeof(char));
  tr.tree1 = (char*)exa_calloc((size_t)tr.treeStringLength, sizeof(char));

  tr.td[0].count = 0;
  tr.td[0].ti    = (traversalInfo *)exa_malloc(sizeof(traversalInfo) * (size_t)tr.mxtips);
  tr.td[0].executeModel = (boolean *)exa_malloc(sizeof(boolean) * (size_t)NUM_BRANCHES);
  tr.td[0].parameterValues = (double *)exa_malloc(sizeof(double) * (size_t)NUM_BRANCHES);

  // tr.fracchange = -1.0;

  tr.constraintVector = (int *)exa_malloc((2 * (size_t)tr.mxtips) * sizeof(int));

  p0 = (nodeptr)exa_malloc((tips + 3 * inter) * sizeof(node));
  assert(p0);

  tr.nodeBaseAddress = p0;


  tr.nodep = (nodeptr *) exa_malloc((2* (size_t)tr.mxtips) * sizeof(nodeptr));
  assert(tr.nodep);    

  tr.nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (int i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash = std::hash<int>()(i); /* hast table stuff */
      p->x      =  0;
      p->xBips  = 0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;
      tr.nodep[i] = p;
    }

  for (int i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (int j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    {
	      p->xBips = 1;
	      p->x = 1;
	    }
	  else
	    {
	      p->xBips = 0;
	      p->x =  0;
	    }
	  p->number = i;
	  p->next   = q;
	  p->back   = (node *) NULL;
	  p->hash   = 0;       
	  q = p;
	}
      p->next->next->next = p;
      tr.nodep[i] = p;
    }

  tr.likelihood  = unlikely;
  tr.start       = (node *) NULL;  

  tr.ntips       = 0;
  tr.nextnode    = 0;
}



void TreeInitializer::initializeBranchLengths(tree &tr, nat numPart, nat numTax )
{
  tr.mxtips = numTax; 
  NUM_BRANCHES = numPart;
  tr.zqr = new double[NUM_BRANCHES]; 
  tr.currentZQR = new double[NUM_BRANCHES]; 
  tr.currentLZR = new double[NUM_BRANCHES];
  tr.currentLZQ = new double[NUM_BRANCHES];
  tr.currentLZS = new double[NUM_BRANCHES];
  tr.currentLZI = new double[NUM_BRANCHES];
  tr.lzs = new double[NUM_BRANCHES];
  tr.lzq = new double[NUM_BRANCHES];
  tr.lzr = new double[NUM_BRANCHES];
  tr.lzi = new double[NUM_BRANCHES];
  tr.coreLZ = new double[NUM_BRANCHES]; 
  tr.curvatOK = new boolean[NUM_BRANCHES];  

  tr.partitionSmoothed = new boolean[NUM_BRANCHES];
  tr.partitionConverged = new boolean[NUM_BRANCHES]; 

  for(nat i = 0; i < numTax ; ++i)
    {
      tr.td[0].ti[i].qz  = new double[NUM_BRANCHES]; 
      tr.td[0].ti[i].rz  = new double[NUM_BRANCHES]; 
    }

  for(nat i = 0; i < numTax + 3 * (numTax - 1)  ; ++i )
    {
      auto node = tr.nodeBaseAddress + i; 
      node->z = new double[NUM_BRANCHES]; 
    }
}

void TreeInitializer::initializeWithAlignmentInfo(TreeAln &traln, RunModes flags)
{
  auto &tr= traln.getTrHandle(); 
  
  traln.setMode(flags); 
  if( ( flags & RunModes::PARTITION_DISTRIBUTION)  !=  RunModes::NOTHING)
    tr.manyPartitions = TRUE; 
  if( (flags & RunModes::MEMORY_SEV) != RunModes::NOTHING)
    tr.saveMemory = TRUE; 

  unifiedInitializePartitions(traln ); 
  initializeBranchLengths(tr, traln.getNumberOfPartitions(), traln.getNumberOfTaxa());
}
