#ifndef WEIBULL_PROPOSER_HPP
#define WEIBULL_PROPOSER_HPP

#include "TreeAln.hpp"
#include "BranchLength.hpp"
#include "extensions.hpp"

class Randomness;

///////////////////////////////////////////////////////////////////////////////
//                              WEIBULL PROPOSER                             //
///////////////////////////////////////////////////////////////////////////////
class WeibullProposer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    WeibullProposer(
        double nrOpt = 0,
        double nrd1 = 0,
        double nrd2 = 0,
        double convParameter = 0,
        double nonConvParameter = 0);
    // ________________________________________________________________________
    BranchLength                           proposeBranch(
        BranchPlain        b,
        TreeAln&           traln,
        AbstractParameter* param,
        Randomness&        rand) const;
    // ________________________________________________________________________
    log_double                             getLogProbability(
        double val) const;
    // ________________________________________________________________________
    bool                                   isConverged() const
    {return _nrD2 < 0 && _nrD1 < 1; }
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&         out,
        const WeibullProposer&rhs)
    {
        out << SHOW(rhs._nrOpt) << SHOW(rhs._nrD1) << SHOW(rhs._nrD2) << SHOW(
            rhs._lambda) << SHOW(rhs._k);
        return out;
    }
    // ________________________________________________________________________
    static std::string                     getName()
    {return std::string("Weibull"); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _nrOpt;      // branch length optimum
    double _nrD1;       // first derivative
    double _nrD2;       // second derivative

    double _lambda;     // shape parameter
    double _k;          // scale parameter

    double _convParameter;
    double _nonConvParameter;
};


#endif
