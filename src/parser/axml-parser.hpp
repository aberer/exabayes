/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#include <stdint.h>

typedef  int boolean;


typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct ent *next;
};

typedef struct ent entry;

typedef unsigned int hashNumberType;

typedef unsigned int parsimonyNumber;

/*typedef uint_fast32_t parsimonyNumber;*/



/*
  typedef uint64_t parsimonyNumber;

  #define PCF 16


typedef unsigned char parsimonyNumber;

#define PCF 2
*/

typedef struct
{
  hashNumberType tableSize;
  entry **table;
  hashNumberType entryCount;
}
  hashtable;


struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;
 
typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;


typedef struct
{
  unsigned int  parsimonyScore;
  unsigned int  parsimonyState;
}
  parsimonyVector;


typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;

// #define NUM_BRANCHES 1 
// #define NUM_PROT_MODELS 21

typedef struct
{
  int tipCase;
  int pNumber;
  int qNumber;
  int rNumber;
  double qz[1];
  double rz[1];
} traversalInfo;

typedef struct
{
  traversalInfo *ti;
  int count;
  int functionType;
  boolean traversalHasChanged;
  boolean *executeModel;
  double  *parameterValues;
} traversalData;


struct noderec;

typedef struct epBrData
{
  int    *countThem;
  int    *executeThem;
  unsigned int *parsimonyScores;
  double *branches;
  double *bootstrapBranches;
  double *likelihoods;
  double originalBranchLength;
  char branchLabel[64];
  int leftNodeNumber;
  int rightNodeNumber;
  int *leftScaling;
  int *rightScaling;
  parsimonyVector *leftParsimony;
  parsimonyVector *rightParsimony;
  double branchLengths[1];
  double *left;
  double *right;
  int branchNumber; 
} epaBranchData;

typedef struct
{
  epaBranchData *epa;

  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;








typedef struct
{
  boolean valid;
  int partitions;
  int *partitionList;
}
  linkageData;

typedef struct
{
  int entries;
  linkageData* ld;
}
  linkageList;


typedef  struct noderec
{
 
  branchInfo      *bInf;
  double           z[1];
#ifdef _BAYESIAN 
  double           z_tmp[1];
#endif 
  struct noderec  *next;
  struct noderec  *back;
  hashNumberType   hash;
  int              support;
  int              number;
  char             x;
}
  node, *nodeptr;

typedef struct
  {
    double lh;
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;




typedef  struct
{
  int              numsp;
  int              sites;
  unsigned char             **y;
  unsigned char             *y0; 
  int              *wgt;
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
} cruncheddata;




typedef struct {
  int     states;
  int     maxTipStates;
  int     lower;
  int      upper;
  int     width;
  int     dataType;
  int     protModels;
  int     autoProtModels;
  int     protFreqs;
  int             **expVector;
  double          **xVector;
  size_t           *xSpaceVector;
 
  unsigned char            **yVector;
  char   *partitionName;
  double *sumBuffer;
 
  double *gammaRates;

  double *EIGN;
  double *EV;



  double *EI;



  

  double *left;
  double *right;




  double *frequencies;
  double *tipVector; 
  double *substRates;
  
  
  double *perSiteRates;

  double *wr;
  double *wr2;

  

  unsigned int    *globalScaler;
  double          *globalScalerDouble; 
  int    *wgt;
 
  int    *rateCategory;
  int    *symmetryVector;
  int    *frequencyGrouping;
  boolean nonGTR;
  boolean optimizeBaseFrequencies;
  double alpha;
  

  int gapVectorLength;
  unsigned int *gapVector;
  double *gapColumn;

  int    numberOfCategories;

  parsimonyNumber parsimonyLength ; 
  parsimonyNumber *parsVect;
} pInfo;



typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;




typedef struct {
 
  int state;

  unsigned int vLength;
  
  int rearrangementsMax;
  int rearrangementsMin;
  int thoroughIterations;
  int fastIterations;
  int treeVectorLength;  
  int mintrav;
  int maxtrav;
  int bestTrav;
  int    Thorough;
  int    optimizeRateCategoryInvocations;
  
  double accumulatedTime;

  double startLH; 
  double lh;
  double previousLh;
  double difference;
  double epsilon;
  
  boolean impr;
  boolean cutoff;  
       
  double tr_startLH;
  double tr_endLH;
  double tr_likelihood;
  double tr_bestOfNode;
  
  double tr_lhCutoff;
  double tr_lhAVG;
  double tr_lhDEC;
  int    tr_NumberOfCategories;
  int    tr_itCount;  
  int    tr_doCutoff;

                                                                    
} checkPointState;


typedef struct {
  double EIGN[19] __attribute__ ((aligned (BYTE_ALIGNMENT)));             
  double EV[400] __attribute__ ((aligned (BYTE_ALIGNMENT)));                
  double EI[380] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double substRates[190];        
  double frequencies[20] ;      
  double tipVector[460] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double fracchange[1];
  double left[1600] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double right[1600] __attribute__ ((aligned (BYTE_ALIGNMENT)));
} siteAAModels;

typedef  struct  {
  boolean useGappedImplementation;
  boolean saveMemory;
  
  siteAAModels siteProtModel[2 * (21  - 2)];// NUM_PROT_MODELS

  boolean estimatePerSiteAA;

  int    *resample;

  int numberOfBranches;
  int    numberOfTipsForInsertion;
  int    *inserts;
  int    branchCounter;

 


  


  parsimonyNumber **parsimonyState_A;
  parsimonyNumber **parsimonyState_C;
  parsimonyNumber **parsimonyState_G;
  parsimonyNumber **parsimonyState_T;
  unsigned int *parsimonyScore; 
  int *ti;
  unsigned int compressedWidth;
  
  int numberOfTrees; 

  stringHashtable  *nameHash;

  pInfo            *partitionData;
  pInfo            *initialPartitionData;
  pInfo            *extendedPartitionData;

  int              *dataVector;
  int              *initialDataVector;
  int              *extendedDataVector;

  int              *patternPosition;
  int              *columnPosition;

  char             *secondaryStructureInput;

  boolean          *executeModel;

  double           *perPartitionLH;

  traversalData    td[1];

  int              maxCategories;

  double           *wr;
  double           *wr2;
  
  double           coreLZ[1];
  int              modelNumber;
  int              numBranches;
  int              bootStopCriterion;
  int              consensusType;
  double           wcThreshold;


 
 
 
 
  
 
  branchInfo	   *bInf;

  int              multiStateModel;


  boolean curvatOK[1];
  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */

  
  double           *fracchanges;

  /* model stuff end */

  unsigned char             **yVector;
  int              secondaryStructureModel;
  uint64_t              originalCrunchedLength;
  int              fullSites;
  int              *originalModel;
  int              *originalDataVector;
  int              *originalWeights;
  int              *secondaryStructurePairs;


  double            *partitionContributions;
  double            fracchange;
  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;

  double           startLH;
  double           endLH;
  double           likelihood;
  double          *likelihoods;
 
  node           **nodep;
  nodeptr          nodeBaseAddress;
  node            *start;
  int              mxtips;  
  int              *model;

  int              *constraintVector;
  int              numberOfSecondaryColumns;
  boolean          searchConvergenceCriterion;
  int              ntips;
  int              nextnode;  
  int              NumberOfModels;
  int              parsimonyLength;
  
  int              checkPointCounter;
  int              treeID;  
  boolean          bigCutoff;
  boolean          partitionSmoothed[1];
  boolean          partitionConverged[1];
  boolean          rooted;
  boolean          grouped;
  boolean          constrained;
  boolean          doCutoff;
  boolean          catOnly;
  rawdata         *rdta;
  cruncheddata    *cdta;

  char **nameList;
  char *tree_string;
  char *tree0;
  char *tree1;
  int treeStringLength;
  unsigned int bestParsimony;
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;

  double zqr[1];
  double currentZQR[1];

  double currentLZR[1];
  double currentLZQ[1];
  double currentLZS[1];
  double currentLZI[1];
  double lzs[1];
  double lzq[1];
  double lzr[1];
  double lzi[1];

 
  int mr_thresh;


  unsigned int **bitVectors;

  unsigned int vLength;

  hashtable *h;
  

} tree;


/***************************************************************/

typedef struct {
  int partitionNumber;
  int partitionLength;
} partitionType;

typedef struct
{
  double z[1];
  nodeptr p, q;
  int cp, cq;
}
  connectRELL, *connptrRELL;

typedef  struct
{
  connectRELL     *connect; 
  int             start;
  double          likelihood;
}
  topolRELL;


typedef  struct
{
  int max;
  topolRELL **t;
}
  topolRELL_LIST;


/**************************************************************/



typedef struct conntyp {
    double           z[1];           /* branch length */
    node            *p, *q;       /* parent and child sectors */
    void            *valptr;      /* pointer to value of subtree */
    int              descend;     /* pointer to first connect of child */
    int              sibling;     /* next connect from same parent */
    } connect, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    connect         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */

    } topol;

typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    boolean          improved;
    } bestlist;

typedef  struct {
  int              categories;
  int              model;
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  boolean          initialSet;
  int              mode;
  long             boot;
  long             rapidBoot;
  boolean          bootstrapBranchLengths;
  boolean          restart;
  boolean          useWeightFile;
  boolean          useMultipleModel;
  boolean          constraint;
  boolean          grouping;
  boolean          randomStartingTree;
  boolean          useInvariant;
  int            protEmpiricalFreqs;
  int            proteinMatrix;
  int            checkpoints;
  int            startingTreeOnly;
  int            multipleRuns;
  long           parsimonySeed;
  boolean        perGeneBranchLengths;
  boolean        likelihoodTest;
  boolean        permuteTreeoptimize;
  boolean        allInOne;
  boolean        generateBS;
  boolean        bootStopping;
  boolean        useExcludeFile;
  boolean        userProteinModel;
  boolean        computeELW;
  boolean        computeDistance;
  boolean        thoroughInsertion;
  boolean        compressPatterns;
  boolean        useSecondaryStructure; 
  double         likelihoodEpsilon;
  double         gapyness;
  int            similarityFilterMode;
  double        *externalAAMatrix; 
  boolean       readTaxaOnly;
  int           meshSearch;  
  boolean       veryFast;
  boolean       useCheckpoint;
  boolean       leaveDropMode;
  int           slidingWindowSize;
  boolean       writeBinaryFile;
  boolean       readBinaryFile;
#ifdef _BAYESIAN 
  boolean       bayesian;
#endif
} analdef;

typedef struct 
{
  int leftLength;
  int rightLength;
  int eignLength;
  int evLength;
  int eiLength;
  int substRatesLength;
  int frequenciesLength;
  int tipVectorLength;
  int symmetryVectorLength;
  int frequencyGroupingLength;

  boolean nonGTR;
  boolean optimizeBaseFrequencies;

  int undetermined;

  const char *inverseMeaning;

  int states;

  boolean smoothFrequencies;

  const unsigned  int *bitVector;

} partitionLengths;

/****************************** FUNCTIONS ****************************************************/



extern void computePlacementBias(tree *tr, analdef *adef);

extern int lookupWord(char *s, stringHashtable *h);

// extern void getDataTypeString(tree *tr, int model, char typeOfData[1024]);

extern unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength);
extern int countTips(nodeptr p, int numsp);
extern entry *initEntry(void);
extern void computeRogueTaxa(tree *tr, char* treeSetFileName, analdef *adef);
extern unsigned int precomputed16_bitcount(unsigned int n);







extern size_t discreteRateCategories(int rateHetModel);

// extern partitionLengths * getPartitionLengths(pInfo *p);
// extern boolean getSmoothFreqs(int dataType);
// extern const unsigned int *getBitVector(int dataType);
// extern int getUndetermined(int dataType);
// extern int getStates(int dataType);
extern char getInverseMeaning(int dataType, unsigned char state);
extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );

extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
// extern boolean whitechar ( int ch );
extern void errorExit ( int e );
extern void printResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBootstrapResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBipartitionResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printLog ( tree *tr, analdef *adef, boolean finalPrint );
extern void printStartingTree ( tree *tr, analdef *adef, boolean finalPrint );
extern void writeInfoFile ( analdef *adef, tree *tr, double t );
extern int main ( int argc, char *argv[] );
extern void calcBipartitions ( tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName );
extern void initReversibleGTR (tree *tr, int model);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void makeGammaCats (double alpha, double *gammaRates, int K);
extern void initModel ( tree *tr, rawdata *rdta, cruncheddata *cdta, analdef *adef );
extern void doAllInOne ( tree *tr, analdef *adef );

extern void classifyML(tree *tr, analdef *adef);
extern void doBootstrap ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void doInference ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void resetBranches ( tree *tr );
extern void modOpt ( tree *tr, analdef *adef , double likelihoodEpsilon);


extern void parsePartitions ( analdef *adef, rawdata *rdta, tree *tr);
extern void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( tree *tr, analdef *adef );
extern void thoroughOptimization ( tree *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( tree *tr, int mintrav, int maxtrav);

extern int checker ( tree *tr, nodeptr p );
extern int randomInt ( int n );
extern void makePermutation ( int *perm, int n, analdef *adef );
extern boolean tipHomogeneityChecker ( tree *tr, nodeptr p, int grouping );
extern void makeRandomTree ( tree *tr, analdef *adef );
extern void nodeRectifier ( tree *tr );
extern void makeParsimonyTreeThorough(tree *tr, analdef *adef);
extern void makeParsimonyTree ( tree *tr, analdef *adef );
extern void makeParsimonyTreeFastDNA(tree *tr, analdef *adef);
extern void makeParsimonyTreeIncomplete ( tree *tr, analdef *adef );
extern void makeParsimonyInsertions(tree *tr, nodeptr startNodeQ, nodeptr startNodeR);



extern FILE *myfopen(const char *path, const char *mode);


extern boolean initrav ( tree *tr, nodeptr p );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern boolean update ( tree *tr, nodeptr p );
extern boolean smooth ( tree *tr, nodeptr p );
extern boolean smoothTree ( tree *tr, int maxtimes );
extern boolean localSmooth ( tree *tr, nodeptr p, int maxtimes );
extern boolean localSmoothMulti(tree *tr, nodeptr p, int maxtimes, int model);
extern void initInfoList ( int n );
extern void freeInfoList ( void );
extern void insertInfoList ( nodeptr node, double likelihood );
extern boolean smoothRegion ( tree *tr, nodeptr p, int region );
extern boolean regionalSmooth ( tree *tr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( tree *tr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( tree *tr, nodeptr p );
extern boolean insertBIG ( tree *tr, nodeptr p, nodeptr q, int numBranches);
extern boolean insertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern boolean testInsertBIG ( tree *tr, nodeptr p, nodeptr q );
extern void addTraverseBIG ( tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( tree *tr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern double treeOptimizeRapid ( tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt);
extern boolean testInsertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( tree *tr );
extern int determineRearrangementSetting ( tree *tr, analdef *adef, bestlist *bestT, bestlist *bt );
extern void computeBIGRAPID ( tree *tr, analdef *adef, boolean estimateModel);
extern boolean treeEvaluate ( tree *tr, double smoothFactor );
extern boolean treeEvaluatePartition ( tree *tr, double smoothFactor, int model );

extern void meshTreeSearch(tree *tr, analdef *adef, int thorough);

extern void initTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, tree *tr, int index );

extern int  saveBestTree (bestlist *bt, tree *tr);
extern int  recallBestTree (bestlist *bt, int rank, tree *tr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern boolean freeBestTree ( bestlist *bt );


extern char *Tree2String ( char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, 
			   boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport);
extern void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission);



extern int treeReadLen (FILE *fp, tree *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly);
extern void treeReadTopologyString(char *treeString, tree *tr);
extern boolean treeReadLenMULT ( FILE *fp, tree *tr, analdef *adef );

extern void getStartingTree ( tree *tr);
extern double treeLength(tree *tr, int model);

extern void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef);
extern boolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength);
extern void computeConsensusOnly(tree *tr, char* treeSetFileName, analdef *adef);
extern double evaluatePartialGeneric (tree *, int i, double ki, int _model);
extern void evaluateGeneric (tree *tr, nodeptr p, boolean fullTraversal);
extern void newviewGeneric (tree *tr, nodeptr p, boolean masked);
extern void newviewGenericMulti (tree *tr, nodeptr p, int model);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern double evaluateGenericVector (tree *tr, nodeptr p);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern boolean isTip(int number, int maxTips);
extern void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches, boolean partialTraversal);



extern void   newviewIterative(tree *tr, int startIndex);

extern void evaluateIterative(tree *);

extern void *malloc_aligned( size_t size);

extern void storeExecuteMaskInTraversalDescriptor(tree *tr);
extern void storeValuesInTraversalDescriptor(tree *tr, double *value);
extern void myBinFwrite(const void *ptr, size_t size, size_t nmemb);
extern void myBinFread(void *ptr, size_t size, size_t nmemb);



extern void makenewzIterative(tree *);
extern void execCore(tree *, volatile double *dlnLdlz, volatile double *d2lnLdlz2);



extern void determineFullTraversal(nodeptr p, tree *tr);
/*extern void optRateCat(tree *, int i, double lower_spacing, double upper_spacing, double *lhs);*/

extern unsigned int evaluateParsimonyIterative(tree *);
extern void newviewParsimonyIterative(tree *);

extern unsigned int evaluateParsimonyIterativeFast(tree *);
extern void newviewParsimonyIterativeFast(tree *);

extern unsigned int evaluatePerSiteParsimony(tree *tr, nodeptr p, unsigned int *siteParsimony);
extern void initravParsimonyNormal(tree *tr, nodeptr p);

extern double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model);
extern void evaluateGenericVectorIterative(tree *, int startIndex, int endIndex);
extern void categorizeIterative(tree *, int startIndex, int endIndex);

extern void fixModelIndices(tree *tr, int endsite, boolean fixRates);
extern void calculateModelOffsets(tree *tr);
extern void gammaToCat(tree *tr);
extern void catToGamma(tree *tr, analdef *adef);
extern void handleExcludeFile(tree *tr, analdef *adef, rawdata *rdta);

extern nodeptr findAnyTip(nodeptr p, int numsp);

extern void parseProteinModel(analdef *adef);



extern void computeNextReplicate(tree *tr, long *seed, int *originalRateCategories, int *originalInvariant, boolean isRapid, boolean fixRates);
/*extern void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);*/

extern void putWAG(double *ext_initialRates);

extern void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant);
extern void parseSecondaryStructure(tree *tr, analdef *adef, int sites);
extern void printPartitions(tree *tr);
extern void compareBips(tree *tr, char *bootStrapFileName, analdef *adef);
extern void computeRF(tree *tr, char *bootStrapFileName, analdef *adef);


extern  unsigned int **initBitVector(tree *tr, unsigned int *vectorLength);
extern hashtable *copyHashTable(hashtable *src, unsigned int vectorLength);
extern hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(hashtable *h, int state);
extern double convergenceCriterion(hashtable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
extern void freeHashTable(hashtable *h);



extern void printBothOpen(const char* format, ... );
extern void printBothOpenMPI(const char* format, ... );
extern void initRateMatrix(tree *tr);

extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf,
				    int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF);

extern int getIncrement(tree *tr, int model);

extern void fastSearch(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta);
extern void shSupports(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta);

extern FILE *getNumberOfTrees(tree *tr, char *fileName, analdef *adef);

extern void writeBinaryModel(tree *tr);
extern void readBinaryModel(tree *tr);
extern void treeEvaluateRandom (tree *tr, double smoothFactor);
extern void treeEvaluateProgressive(tree *tr);

extern void testGapped(tree *tr);

extern boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen);
extern boolean compatible(entry* e1, entry* e2, unsigned int bvlen);



extern int *permutationSH(tree *tr, int nBootstrap, long _randomSeed);

extern void updatePerSiteRates(tree *tr, boolean scaleRates);

extern void restart(tree *tr);







