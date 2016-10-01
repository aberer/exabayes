#ifndef _PROPOSAL_FUNCTION_H
#define _PROPOSAL_FUNCTION_H

#include "Randomness.hpp"
#include "Density.hpp"
#include "Chain.hpp"
#include "AbstractProposal.hpp"

#include <vector>
#include <algorithm>

///////////////////////////////////////////////////////////////////////////////
//                             ABSTRACT PROPOSER                             //
///////////////////////////////////////////////////////////////////////////////
class AbstractProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC TYPES                             //
    ///////////////////////////////////////////////////////////////////////////
public:
    using UPtr = std::unique_ptr<AbstractProposer>;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    AbstractProposer(
        bool   tune,
        bool   tuneup,
        double minVal,
        double maxVal);
    // ________________________________________________________________________
    virtual ~AbstractProposer(){}
    // ________________________________________________________________________
    virtual std::vector<double>                    proposeValues(
        std::vector<double>oldValues,
        double             parameter,
        Randomness&        rand,
        log_double&        hastings) = 0;
    // ________________________________________________________________________
    bool                                           isTune() const
    {return _tune;}
    // ________________________________________________________________________
    bool                                           isTuneup() const
    {return _tuneup;}
    // ________________________________________________________________________
    void                                           correctAbsoluteRates(
        std::vector<double>&values) const;
    // ________________________________________________________________________
    virtual AbstractProposer*                      clone() const = 0;

    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    bool _tune;
    bool _tuneup;
    double _minVal;
    double _maxVal;
};


#endif

