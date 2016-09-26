#ifndef _PARTITION_LIKELIHOOD
#define _PARTITION_LIKELIHOOD

#include "TreeAln.hpp"

#include <vector>

///////////////////////////////////////////////////////////////////////////////
//                            PARTITION LIKELIHOOD                           //
///////////////////////////////////////////////////////////////////////////////
class PartitionLikelihood
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    nat                  model;
    std::vector<double*> cachedArrays;
    std::vector<size_t>  lengths;
    std::vector<nat>     scaler;
    std::vector<bool>    isCached;
    std::vector<nat>     gapVector;
    std::vector<double>  gapColumn;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    PartitionLikelihood(
        const TreeAln& traln,
        nat            modelArg,
        bool           useSEV)
        : model(modelArg)
        , cachedArrays(traln.getNumberOfInnerNodes(false))
        ,                                              // we are talking about
                                                       // the unrooted pll tree
        lengths(traln.getNumberOfInnerNodes(false), 0)
        , scaler(2 * traln.getNumberOfTaxa(), 0)
        , isCached(traln.getNumberOfInnerNodes(false), false)
        , gapVector{}
        , gapColumn{}
    {
        if (useSEV)
        {
// #if 0
            auto& partition = traln.getPartition(model);
            gapVector = std::vector<nat>(partition.getGapVectorLength() * 2
                                         * traln.getNumberOfTaxa(), 0);
            gapColumn = std::vector<double>(traln.getNumberOfTaxa()
                                            * partition.getStates()
                                            * Partition::maxCategories, 0);
            // very MEH

// #else
//  // TODO gap stuff
//  assert(0);
// #endif
        }
    }
};


#endif
