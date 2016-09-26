#ifndef _OPTIMIZED_PARAMETER_HPP
#define _OPTIMIZED_PARAMETER_HPP

#include "AbstractParameter.hpp"

#include "BranchPlain.hpp"
#include "DistributionProposer.hpp"

class TreeAln;


///////////////////////////////////////////////////////////////////////////////
//                            OPTIMIZED PARAMETER                            //
///////////////////////////////////////////////////////////////////////////////
class OptimizedParameter
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
    // TODO: really?
public:
    static const double zMin;
    static const double zMax;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    OptimizedParameter(
        TreeAln&           traln,
        const BranchPlain& branch,
        AbstractParameter* param,
        int                maxIter);
    // ________________________________________________________________________
    template<class T>
    auto                                   getProposerDistribution(
        TreeAln&traln,
        double  convParameter,
        double  nonConvParameter) const
        ->DistributionProposer<T>;
    // ________________________________________________________________________
    void                                   applyToMask(
        std::vector<bool>&mask)  const;
    // ________________________________________________________________________
    void                                   applyValues(
        double*values) const;
    // ________________________________________________________________________
    BranchLength                           getOptimizedBranch() const;
    // ________________________________________________________________________
    AbstractParameter*                     getParam() const
    {return _param; }
    // ________________________________________________________________________
    bool                                   isCurvatureOk() const
    {return _curvatOK; }
    // ________________________________________________________________________
    bool                                   hasConvergedNew() const;
    // ________________________________________________________________________
    bool                                   hasFinished() const;
    // ________________________________________________________________________
    void                                   resetStep();
    // ________________________________________________________________________
    void                                   changeSide();
    // ________________________________________________________________________
    void                                   decrIter()
    {--_iters; }
    // ________________________________________________________________________
    void                                   shortenBadBranch();
    // ________________________________________________________________________
    void                                   nrStep();
    // ________________________________________________________________________
    void                                   setToTypicalBranch(
        double   typicalAbsLen,
        TreeAln& traln);
    // ________________________________________________________________________
    void                                   extractDerivatives(
        TreeAln&            traln,
        std::vector<double>&dlnLdlz,
        std::vector<double>&d2lnLdlz2);
    // ________________________________________________________________________
    // check if we are done for partition i or if we need to adapt the
    // branch length again
    //
    void                                   doInitStep();

    // ________________________________________________________________________
    double                                 getFirstDerivative() const
    {return _nrD1; }
    // ________________________________________________________________________
    double                                 getSecondDerivative()
    const
    {return _nrD2; }
    // ________________________________________________________________________
    double                                 getOptimum() const
    {return _zCur; }
    // ________________________________________________________________________
    void                                   checkConvergence();
    // ________________________________________________________________________
    BranchPlain                            getBranch() const
    {return _branch; }
    // ________________________________________________________________________
    friend std::ostream&                   operator                    <<(
        std::ostream&             out,
        const OptimizedParameter& rhs)
    {
        out << rhs.getOptimum() << "," << rhs.getFirstDerivative() << ","
            << rhs.getSecondDerivative();
        return out;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    double             _zPrev;
    double             _zCur;
    double             _zStep;
    double             _coreLZ;

    double             _nrD1;
    double             _nrD2;

    int                _iters;

    bool               _curvatOK;
    AbstractParameter* _param;

    bool               _hasConverged;
    BranchPlain        _branch;
};


///////////////////////////////////////////////////////////////////////////////
//                             INLINE DEFINITIONS                            //
///////////////////////////////////////////////////////////////////////////////
template<class T>
auto                                       OptimizedParameter::
    getProposerDistribution(
    TreeAln&traln,
    double  convParameter,
    double  nonConvParameter) const
    ->DistributionProposer<T>
{
    // tout << "init with factor "  << factor << std::endl;
    auto bl =  BranchLength(BranchPlain(1, 2));
    bl.setLength(_zCur);
    auto realLen = bl.toMeanSubstitutions(_param->getMeanSubstitutionRate());
    return DistributionProposer<T>(realLen, _nrD1, _nrD2, convParameter,
                                   nonConvParameter);
}


#endif

