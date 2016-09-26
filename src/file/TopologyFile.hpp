#ifndef TOPOLOGY_FILE
#define TOPOLOGY_FILE

#include "common.h"
#include "TreeAln.hpp"
#include "OutputFile.hpp"

#include <sstream>
#include <string>
#include <iostream>

///////////////////////////////////////////////////////////////////////////////
//                               TOPOLOGY FILE                               //
///////////////////////////////////////////////////////////////////////////////
class TopologyFile : public OutputFile
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    TopologyFile(
        std::string workdir,
        std::string runname,
        nat         runid,
        nat         couplingId,
        nat         paramNum,
        bool        hasManyTopoloFiles);
    // ________________________________________________________________________
    void                              initialize(
        const TreeAln& traln,
        nat            someId,
        bool           isDryRun);
    // ________________________________________________________________________
    void                              sample(
        const TreeAln&     traln,
        uint64_t           gen,
        AbstractParameter* blParams);
    // ________________________________________________________________________
    void                              regenerate(
        std::string workdir,
        std::string prevId,
        uint64_t    gen);
    // ________________________________________________________________________
    void                              verifyNonExistance();

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::streamoff                    getPosBeforeEnd() const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    nat  runid;
    nat  couplingId;
    nat  paramNum;
    bool hasManyTopoloFiles;
};


#endif
