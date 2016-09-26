#ifndef CONSENSUSTREE_HPP
#define CONSENSUSTREE_HPP

#include "BipartitionExtractor.hpp"
#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////
//                               CONSENSUS TREE                              //
///////////////////////////////////////////////////////////////////////////////
class ConsensusTree
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ConsensusTree(
        std::vector<std::string>files,
        double                  burnin,
        double                  threshold,
        bool                    isMRE);
    // ________________________________________________________________________
    std::string                                 getConsensusTreeString(
        bool printNames) const;
    // ________________________________________________________________________
    std::vector<Bipartition>                    getRefinedConsensusTree(
        const std::vector<Bipartition>&consensusBips,
        const std::vector<Bipartition>&minorityBips) const;
    // ________________________________________________________________________
    std::string                                 getTreeHeader() const;
    // ________________________________________________________________________
    std::string                                 getType() const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    BipartitionExtractor bipEx;
    double               _threshold;
    bool                 _isMRE;
};


#endif
