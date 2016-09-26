
/**
 * @file ConfigReader.hpp
 *
 * @brief specializes a nxsreader for parsing of an exabayes block.
 */

#ifndef _NCL_CONFIG_READER
#define _NCL_CONFIG_READER

#include <ncl/ncl.h>
#include <vector>

#include "PriorBelief.hpp"
#include "AbstractProposal.hpp"

#include "BlockParams.hpp"
#include "BlockPrior.hpp"

///////////////////////////////////////////////////////////////////////////////
//                               CONFIG READER                               //
///////////////////////////////////////////////////////////////////////////////
class ConfigReader : public NxsReader
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ConfigReader()
        : NxsReader(){SetWarningOutputLevel(SUPPRESS_WARNINGS_LEVEL); }
    // ________________________________________________________________________
    virtual void                    ExitingBlock(
        NxsString blockName){}
    // ________________________________________________________________________
    virtual void                    ExecuteStopping(){}
    // ________________________________________________________________________
    virtual void                    ExecuteStarting(){}
};


#endif
