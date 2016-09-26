#ifndef BIPARTITION_EXTRACTOR_HPP
#define BIPARTITION_EXTRACTOR_HPP

#include "TreeProcessor.hpp"
#include "BipartitionHash.hpp"

#include <unordered_map>

///////////////////////////////////////////////////////////////////////////////
//                           BIPARTITION EXTRACTOR                           //
///////////////////////////////////////////////////////////////////////////////
class BipartitionExtractor : public TreeProcessor
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BipartitionExtractor(
        std::vector<std::string>files,
        bool                    extractToOneHash,
        bool                    expensiveCheck);
    // ________________________________________________________________________
    BipartitionExtractor(
        BipartitionExtractor&& rhs) = delete;
    // ________________________________________________________________________
    BipartitionExtractor&
                        operator=(
        BipartitionExtractor rhs) = delete;
    // ________________________________________________________________________
    virtual ~BipartitionExtractor(){}
    // ________________________________________________________________________
    template<bool readBL>
    void
                        extractBips(
        nat burnin);
    // ________________________________________________________________________
    void
                        printBipartitions(
        std::string id) const;
    // ________________________________________________________________________
    void
                        printBipartitionStatistics(
        std::string id) const;
    // ________________________________________________________________________
    void
                        printFileNames(
        std::string id) const;
    // ________________________________________________________________________
    void
                        printBranchLengths(
        std::string id) const;
    // ________________________________________________________________________
    const std::vector<BipartitionHash>&
                        getBipartitionHashes()
    const {return _bipHashes; }
    // ________________________________________________________________________
    std::string
                        bipartitionsToTreeString(
        std::vector<Bipartition>bips,
        bool                    printSupport,
        bool                    printBranchLengths,
        bool                    phylipStyle) const;
    // ________________________________________________________________________
    nat
                        getNumTreesInFile(
        std::string file) const;
    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void
                        extractUniqueBipartitions();
    // ________________________________________________________________________
    void
                        buildTreeRecursive(
        nat                                   currentId,
        const std::vector<std::vector<nat> >& directSubBips,
        const std::vector<Bipartition>&       bips,
        std::stringstream&                    result,
        bool                                  printSupport,
        bool                                  printBranchLengths,
        bool                                  phylipStyle) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<BipartitionHash>         _bipHashes;
    std::unordered_map<Bipartition, nat> _uniqueBips;
    bool                                 _extractToOneHash;
};


#endif
