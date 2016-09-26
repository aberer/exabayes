#ifndef _DIAGNOSTICS_FILE_H
#define _DIAGNOSTICS_FILE_H

#include "CoupledChains.hpp"
#include "OutputFile.hpp"
#include "GlobalVariables.hpp"

#include <vector>

///////////////////////////////////////////////////////////////////////////////
//                              DIAGNOSTICS FILE                             //
///////////////////////////////////////////////////////////////////////////////
class DiagnosticsFile : public OutputFile
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    DiagnosticsFile()
        : _names{}
        , _initialized(false)
    {}
    // ________________________________________________________________________
    void                           initialize(
        std::string                      workdir,
        std::string                      name,
        const std::vector<CoupledChains>&runs);
    // ________________________________________________________________________
    void                           regenerate(
        std::string workdir,
        std::string nowId,
        std::string prevId,
        uint64_t    gen);
    // ________________________________________________________________________
    void                           printDiagnostics(
        uint64_t                         gen,
        double                           asdsf,
        const std::vector<CoupledChains>&runs);
    // ________________________________________________________________________
    bool                           isInitialized() const
    {return _initialized; }

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    std::string                    createName(
        std::string runname,
        std::string workdir);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::vector<std::string> _names;
    bool                     _initialized;
};


#endif
