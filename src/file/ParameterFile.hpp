#ifndef _PARAMETER_FILE
#define _PARAMETER_FILE

#include "ParameterList.hpp"
#include "AbstractParameter.hpp"
#include "OutputFile.hpp"

#include <sstream>
#include <string>
#include <iostream>
#include <limits>
#include <vector>

/**
 * @notice opening and closing the stream all the time may appear
 * weird. Since properly the copy/move constructor for streams is not
 * implemented properly in g++, it appears as a good choice to do it
 * as implemented below.
 */

///////////////////////////////////////////////////////////////////////////////
//                               PARAMETER FILE                              //
///////////////////////////////////////////////////////////////////////////////
class ParameterFile : public OutputFile
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ParameterFile(
        std::string workdir,
        std::string runname,
        nat         runid);
    // ________________________________________________________________________
    void                    initialize(
        const TreeAln&       traln,
        const ParameterList& parameters,
        nat                  someId,
        bool                 isDryRun);
    // ________________________________________________________________________
    void                    sample(
        const TreeAln&       traln,
        const ParameterList& parameters,
        uint64_t             gen,
        log_double           lnPr);
    // ________________________________________________________________________
    void                    regenerate(
        std::string workdir,
        std::string prevId,
        uint64_t    gen);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    nat runid;
};


#endif
