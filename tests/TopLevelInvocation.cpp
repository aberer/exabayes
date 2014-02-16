#include <gtest/gtest.h>
#include <sstream>

#include <sys/sendfile.h>  // sendfile
#include <fcntl.h>         // open
#include <unistd.h>        // close
#include <sys/stat.h>      // fstat
#include <sys/types.h>     // fstat


#ifdef _WITH_MPI
#include <mpi.h>
#endif

int NUM_BRANCHES; 

#define _INCLUDE_DEFINITIONS
#include "GlobalVariables.hpp"
#undef _INCLUDE_DEFINITIONS

#include  "../src/SampleMaster.hpp"


/** 
    Tests if various command line parameters result in the same known
    likelihood.
*/ 

// that was a tough one =(
extern int optind; 


class TopLevel  : public testing::Test
{
public: 
  void SetUp()
  {
    globals.logFile = "/dev/null";   
    globals.logStream =  new std::ofstream{"/dev/null"}; 
    globals.teeOut =  new teestream(*globals.logStream, *globals.logStream);
  }
}; 


static std::vector<string> convertLine(std::string line)
{
  auto &&iss = std::istringstream{line}; 
  auto result = std::vector<std::string>{}; 
  
  while(iss)
    {
      auto str = std::string{}; 
      iss >> str; 
      if(str.compare("") != 0 )
	result.push_back(str) ; 
    }
  return result; 
}


CommandLine constructCommandLine(const std::string& str)
{  
  auto arr = convertLine(str); 
  auto arr2 = std::vector<char*>{}; 
  for(auto &elem : arr)
    arr2.push_back(strdup(elem.c_str())); 
  auto cl = CommandLine{}; 

  cl.initialize(arr2.size(), arr2.data());
  return  cl; 
}


SampleMaster init( std::string cmdln )
{
  optind = 0; 
  
  auto cl = constructCommandLine(cmdln) ; 
  auto _plPtr = make_shared<ParallelSetup>(); 
#if HAVE_PLL == 0 
  _plPtr->initializeExaml(cl);
#endif

  char** tmp; 
  auto _sm = SampleMaster{}; 
  _sm.setParallelSetup(_plPtr);
  _sm.setCommandLine(cl); 
  _sm.initializeRuns(Randomness(cl.getSeed())); 
  return _sm; 
}


void assertValues(SampleMaster& _sm , double trueBest, double truePp)
{
  auto &c = _sm.getRuns()[0].getChains()[0]; 
  auto best = c.getBestState(); 
  auto pp = c.getLikelihood() + c.getPrior().getLnPrior(); 

  // std::cout << MAX_SCI_PRECISION << "best="  << best << "\tpp=" << pp << std::endl; 

  ASSERT_TRUE( fabs( best - trueBest )  < 1e-6); 
  ASSERT_TRUE( fabs (pp - truePp ) < 1e6-6 ); 
}


double dnaBest = -1.607602001225307e+04; 
double dnaPP = -1.579327993707598e+04; 

double aaBest=-2.950608611893024e+03; 
double aaPP = -2.996126255407966e+03; 


double bestDNA[4] = {  -1.606184140768954e+04 , -1.605852797997213e+04, -1.604100762889591e+04, -1.604441839175527e+04} ; 
void assertMultiValues(SampleMaster& _sm, bool onlyAtMaster )
{
  _sm.synchronize(CommFlag::PRINT_STAT | CommFlag::PROPOSALS | CommFlag::TREE | CommFlag::SWAP);

  int myRank = 0; 
#ifdef _WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank );
#endif
  if(onlyAtMaster && myRank == 0 )
    {
      nat ctr = 0; 
      for(auto &run : _sm.getRuns())
	{
	  for(auto &c : run.getChains() )
	    {
	      ASSERT_TRUE(fabs(c.getBestState()  - bestDNA[ctr] )  < 1e-6);
	      ++ctr; 
	    }
	}
    }
}



static void myCopyFile(std::string src, std::string desti)
{
  int source = open(src.c_str(), O_RDONLY, 0);
  int dest = open(desti.c_str(), O_WRONLY | O_CREAT, 0644);

  struct stat stat_source;
  fstat(source, &stat_source);

  sendfile(dest, source, 0, stat_source.st_size);

  close(source);
  close(dest);
}


static void copyResumeFiles()
{
  auto resumeFolder = std::string{"./test-data/dna/toResume"}; 
  auto files = std::vector<std::string> {"ExaBayes_checkpoint.toResume" ,"ExaBayes_info.toResume", "ExaBayes_parameters.toResume.1", "ExaBayes_topologies.toResume.0.tree.0", "ExaBayes_topologies.toResume.1.tree.0", "ExaBayes_diagnostics.toResume",  "ExaBayes_parameters.toResume.0",  "ExaBayes_prevCheckpointBackup.toResume",  "ExaBayes_topologies.toResume.0.tree.1",  "ExaBayes_topologies.toResume.1.tree.1"}; 
  
  for(auto file : files)
    myCopyFile(std::string{resumeFolder + "/"  + file}, std::string{ "./" + file});
}



////////////////
// THE TESTS  //
////////////////

// DNA parted
TEST_F(TopLevel, dnaPartedBinary)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnaPartedBinary -f " + folder +  "/aln.binary -c " + folder + "/config.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest ,dnaPP);
}


// parsing DNA  
TEST_F(TopLevel, dnaPartedRaw)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init("bla -s 123 -n test-dnaPartedRaw -f " + folder + "/aln.phy -q " + folder + "/aln.model  -c " + folder + "/config.nex");
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest, dnaPP);
}


// AA parted
TEST_F(TopLevel, aaParted)
{  
  auto folder = std::string{"test-data/aa"}; 
  auto _sm = init(  "bla -s 123 -n test-aaParted -f "+ folder+  "/aln.binary -c " + folder + "/config.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  
  assertValues(_sm, aaBest, aaPP);
}


// parsing AA 
TEST_F(TopLevel, aaPartedRaw)
{  
  auto folder = std::string{"test-data/aa"}; 
  auto _sm = init(  "bla -s 123 -n test-aaParted -f " + folder +  "/aln.phy  -q "  + folder + "/aln.model -c " + folder + "/config.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, aaBest, aaPP);
}


// minimal aa (4 taxa)
TEST_F(TopLevel, minimal)
{  
  auto folder = std::string{"test-data/aa-min"}; 
  auto _sm = init(  "bla -s 123 -n test-aamin -f " + folder +  "/aln.phy  -q "  + folder + "/aln.model -c " + folder + "/config.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, -1.838165439052971e+03, -1.882050089193929e+03);
}



// mixed
TEST_F(TopLevel, mixed)
{  
  auto folder = std::string{"test-data/mixed"}; 
  auto _sm = init(  "bla -s 123 -n test-aamin -f " + folder +  "/aln.phy  -q "  + folder + "/aln.model -c " + folder + "/config.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm,-7.797382559009363e+04,-7.694705148165907e+04);
}


// mixed  with SEV vectors
TEST_F(TopLevel, mixedSEV)
{  
  auto folder = std::string{"test-data/mixed"}; 
  auto _sm = init(  "bla -s 123 -n test-aamin-sev -f " + folder +  "/aln.phy  -q "  + folder + "/aln.model -c " + folder + "/config.nex -S "  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm,-7.797382559009363e+04,-7.694705148165907e+04);
}


// test -M 1 
TEST_F(TopLevel, dnamem1)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnamem1 -f " + folder +  "/aln.binary -c " + folder + "/config.nex -M 1"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest ,dnaPP);
}


// test -M 2 
TEST_F(TopLevel, dnamem2)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnamem2 -f " + folder +  "/aln.binary -c " + folder + "/config.nex -M 2"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest ,dnaPP);
}


// test -M 3 
TEST_F(TopLevel, dnamem3)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnamem3 -f " + folder +  "/aln.binary -c " + folder + "/config.nex -M 3" );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest ,dnaPP);
}


// test -M 3 -S 
TEST_F(TopLevel, dnamem3sev)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnamem3sev -f " + folder +  "/aln.binary -c " + folder + "/config.nex -M 3 -S"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm, dnaBest ,dnaPP);
}


// test multiple runs and chains 
TEST_F(TopLevel, dnaMulti)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnaMulti -f " + folder +  "/aln.binary -c " + folder + "/config-multi.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertMultiValues(_sm, false);
}


////////////////////////////
///// NON-MPI specific /////
////////////////////////////
#ifndef _WITH_MPI
// test checkpoint 
TEST_F(TopLevel, dnaMultiCheckpoint)
{  
  auto folder = std::string{"test-data/dna"}; 
  copyResumeFiles();

  auto _sm = init(  "bla -s 123 -n test-dnaMultiCheckpoint -r toResume -f " + folder +  "/aln.binary -c " + folder + "/config-multi.nex"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertMultiValues(_sm, false);
}
#endif


///////////////////
// MPI specific  //
///////////////////

// some tests that only make sense with MPI 
#ifdef _WITH_MPI

// resume from a checkpoint in parallel    
TEST_F(TopLevel, dnaMultiCheckpoint)
{  
  auto folder = std::string{"test-data/dna"}; 
  copyResumeFiles();

  auto _sm = init(  "bla -n test-dnaMultiCheckpoint -r toResume -f " + folder +  "/aln.binary -c " + folder + "/config-multi.nex -C 2 "  );
  _sm.run();
  _sm.deleteMyFiles();
  assertMultiValues(_sm, false);
}



// test multiple runs and chains, 2 chains parallel 
TEST_F(TopLevel, dnaMultiParaChain)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnaMultiParaChain -f " + folder +  "/aln.binary -c " + folder + "/config-multi.nex -C 2"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertMultiValues(_sm, false);
}


// test multiple runs and chains, 2 runs parallel 
TEST_F(TopLevel, dnaMultiParaRun)
{  
  auto folder = std::string{"test-data/dna"}; 
  auto _sm = init(  "bla -s 123 -n test-dnamultipararun -f " + folder +  "/aln.binary -c " + folder + "/config-multi.nex -R 2"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertMultiValues(_sm, true);
}


// mixed
TEST_F(TopLevel, mixedWithQ)
{  
  auto folder = std::string{"test-data/mixed"}; 
  auto _sm = init(  "bla -s 123 -n test-aamin-q -f " + folder +  "/aln.phy  -q "  + folder + "/aln.model -c " + folder + "/config.nex -Q"  );
  _sm.run();
  _sm.deleteMyFiles();
  assertValues(_sm,-7.797382559009363e+04,-7.694705148165907e+04);
}

#endif
