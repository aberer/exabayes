#ifndef _INSERTION_RESULT_HPP
#define  _INSERTION_RESULT_HPP

#include "OptimizedParameter.hpp"

///////////////////////////////////////////////////////////////////////////////
//                              INSERTION RESULT                             //
///////////////////////////////////////////////////////////////////////////////
class InsertionResult
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    InsertionResult(
        BranchPlain branch = BranchPlain(0, 0),
        log_double prob = log_double::fromAbs(1.),
        std::vector<OptimizedParameter>optParams =
            std::vector<OptimizedParameter>{})
        : _branch{branch},
    _lprob{prob},
    _optParams{optParams}
    {}
    // ________________________________________________________________________
    BranchPlain                                        getBranch() const
    {return _branch;}
    // ________________________________________________________________________
    log_double                                         getInsertionProb() const
    {return _lprob;}
    // ________________________________________________________________________
    void                                               setInsertionProb(
        log_double d){_lprob = d;}
    // ________________________________________________________________________
    std::vector<OptimizedParameter>                    getOptParams() const
    {return _optParams;}
    // ________________________________________________________________________
    friend std::ostream&                               operator<<(
        std::ostream&          out,
        const InsertionResult& rhs)
    {
        out << rhs._branch << "\t" << rhs._lprob << "\t"
            << rhs._optParams.size();
        return out;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    BranchPlain _branch;
    log_double _lprob;
    std::vector<OptimizedParameter>_optParams;
};


#endif
