#ifndef _COMMANDLINE_H
#define _COMMANDLINE_H

#include "Randomness.hpp"
#include "MemoryMode.hpp"

#include <memory>
#include <iosfwd>
#include <string>

///////////////////////////////////////////////////////////////////////////////
//                               COMMAND LINE                                //
///////////////////////////////////////////////////////////////////////////////
class CommandLine
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    CommandLine();
    // ________________________________________________________________________
    void                           initialize(
        int   argc,
        char**argv);
    // ________________________________________________________________________
    randCtr_t                      getSeed() const;
    // ________________________________________________________________________
    std::string                    getConfigFileName() const
    {return configFileName; }
    // ________________________________________________________________________
    std::string                    getAlnFileName() const
    {return alnFileName; }
    // ________________________________________________________________________
    std::string                    getRunid() const
    {return runid; }
    // ________________________________________________________________________
    std::string                    getTreeFile() const
    {return treeFile; }
    // ________________________________________________________________________
    int                            getNumRunParallel() const
    {return runNumParallel; }
    // ________________________________________________________________________
    std::string                    getWorkdir() const
    {return workDir; }
    // ________________________________________________________________________
    void                           printVersion(
        std::ostream&out);
    // ________________________________________________________________________
    nat                            getNumChainsParallel() const
    {return chainNumParallel; }
    // ________________________________________________________________________
    std::string                    getCheckpointId() const
    {return checkpointId; }
    // ________________________________________________________________________
    void                           parseAlternative(
        int  argc,
        char*argv[]);
    // ________________________________________________________________________
    bool                           isSaveMemorySEV() const
    {return saveMemorySEV; }
    // ________________________________________________________________________
    std::string                    getCommandLineString() const;
    // ________________________________________________________________________
    MemoryMode                     getMemoryMode() const
    {return memoryMode; }
    // ________________________________________________________________________
    bool                           isDryRun() const
    {return dryRun; }
    // ________________________________________________________________________
    bool                           onlyPrintHelp() const
    {return _onlyPrintHelp; }
    // ________________________________________________________________________
    bool                           onlyPrintVersion() const
    {return _onlyPrintVersion; }
    // ________________________________________________________________________
    bool                           alnFileIsBinary() const;
    // ________________________________________________________________________
    bool                           hasThreadPinning() const
    {return _hasThreadPinning; }
    // ________________________________________________________________________
    int                            getNumThreads() const
    {return _totalThreads; }
    // ________________________________________________________________________
    std::string                    getSingleModel() const
    {return singleModel; }
    // ________________________________________________________________________
    std::string                    getModelFile() const
    {return modelFile; }
    // ________________________________________________________________________
    bool                           isQuiet() const
    {return quiet; }
    // ________________________________________________________________________
    int                            getReaderStride() const
    {return readerStride; }
    // ________________________________________________________________________
    void                           printHelp(
        std::ostream& out);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    void                           assertFileExists(
        std::string filename);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    randCtr_t   seed;
    std::string configFileName;
    std::string alnFileName;
    std::string runid;
    std::string treeFile;
    std::string workDir;
    int         runNumParallel;
    int         chainNumParallel;
    std::string checkpointId;
    MemoryMode  memoryMode;
    bool        saveMemorySEV;
    bool        dryRun;
    std::string modelFile;
    std::string singleModel;
    bool        quiet;
    int         readerStride;
    std::string _cmdString;
    int         _totalThreads;
    bool        _hasThreadPinning;
    bool        _onlyPrintHelp;
    bool        _onlyPrintVersion;
};


#endif
