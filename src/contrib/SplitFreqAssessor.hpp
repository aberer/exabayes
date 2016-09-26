/**
 *  @file SplitFreqAssessor.hpp
 *
 *  @brief Calculates the asdsf.
 *
 *  @notice This file is full of hacks, to get this somehow going.
 */

#ifndef _AVGSPLITFREQASSESSOR_H
#define _AVGSPLITFREQASSESSOR_H

#include "TreeProcessor.hpp"
#include "BipartitionHash.hpp"

///////////////////////////////////////////////////////////////////////////////
//                            SPLIT FREQ ASSESSOR                            //
///////////////////////////////////////////////////////////////////////////////
class SplitFreqAssessor : public TreeProcessor
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    SplitFreqAssessor(
        std::vector<std::string>fileNames,
        bool                    expensiveCheck);
    // ________________________________________________________________________
    SplitFreqAssessor(
        const SplitFreqAssessor&rhs) = default;
    // ________________________________________________________________________
    SplitFreqAssessor(
        SplitFreqAssessor&&rhs) = default;
    // ________________________________________________________________________
    SplitFreqAssessor&                           operator=(
        const SplitFreqAssessor&rhs)  = default;
    // ________________________________________________________________________
    SplitFreqAssessor&                           operator=(
        SplitFreqAssessor&&rhs) = default;
    // ________________________________________________________________________
    virtual ~SplitFreqAssessor(){}
    /**
     *  @brief use the new bipartition hash for extracting bipartitions
     */
    void                                         extractBips(
        const std::vector<nat>&start,
        const std::vector<nat>&end);
    // ________________________________________________________________________
    /**
     *  @brief gets the minimum number of trees present in all of the files
     */
    int                                          getMinNumTrees();
    // ________________________________________________________________________
    std::pair<double, double>                    computeAsdsfNew(
        double ignoreFreq);
    // ________________________________________________________________________
    nat                                          getNumTreeAvailable(
        std::string filename);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<BipartitionHash>         newBipHashes;
    std::unordered_map<std::string, nat> file2numTree;
};


#endif
