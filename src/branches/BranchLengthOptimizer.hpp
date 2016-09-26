#ifndef _BRANCH_LENGTH_OPT_HPP
#define  _BRANCH_LENGTH_OPT_HPP


#include "OptimizedParameter.hpp"

#include "AbstractParameter.hpp"
#include "Communicator.hpp"

#include <vector>
#include <functional>

class TreeAln;

///////////////////////////////////////////////////////////////////////////////
//                          BRANCH LENGTH OPTIMIZER                          //
///////////////////////////////////////////////////////////////////////////////
class BranchLengthOptimizer
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BranchLengthOptimizer(
        TreeAln&                              traln,
        const BranchPlain&                    branch,
        int                                   maxIter,
        Communicator&                         comm,
        const std::vector<AbstractParameter*>&blParams);
    // ________________________________________________________________________
    void
                        applyToTraversalDescriptor(
        std::vector<bool>&execModel,
        TreeAln&          traln) const;
    // ________________________________________________________________________
    bool
                        hasConvergedAll()  const;
    // ________________________________________________________________________
    /**
     *  @brief this function is essentially makenewzGeneric
     */
    void
                        optimizeBranches(
        TreeAln& traln);
    // ________________________________________________________________________
    std::vector<OptimizedParameter>
                        getOptimizedParameters()
    const {return _optParams; }

    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::reference_wrapper<Communicator>  _comm;
    const std::vector<AbstractParameter*> _blParams;
    std::vector<OptimizedParameter>       _optParams;
    BranchPlain                           _branch;
    std::vector<BranchLength>             _origBranches;
};


#endif
