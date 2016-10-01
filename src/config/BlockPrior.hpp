#ifndef _BLOCK_PRIOR_H
#define _BLOCK_PRIOR_H

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "GlobalVariables.hpp"

#include "ExaBlock.hpp"

#include "AbstractPrior.hpp"

#include "Category.hpp"

// if the set is empty, then we have a general "fall-back" prior
using multiMapCategory2TuplePartitionsPrior =
        std::unordered_multimap<Category,
                                std::tuple<std::unordered_set<nat>,
                                           AbstractPrior::UPtr> >;


///////////////////////////////////////////////////////////////////////////////
//                                BLOCK PRIOR                                //
///////////////////////////////////////////////////////////////////////////////
class BlockPrior : public ExaBlock
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit BlockPrior(
        size_t numPart)
        : _parsedPriors{}
        , _numPart(numPart)
    {NCL_BLOCKTYPE_ATTR_NAME = "PRIORS";}
    // ________________________________________________________________________
    void                                                           verify()
    const;
    // ________________________________________________________________________
    virtual void                                                   Read(
        NxsToken&token);
    // ________________________________________________________________________
    const multiMapCategory2TuplePartitionsPrior&                   getPriors()
    const {return _parsedPriors;}

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    AbstractPrior::UPtr                                            parsePrior(
        NxsToken&token);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    multiMapCategory2TuplePartitionsPrior
        _parsedPriors;
    size_t _numPart;
};


#endif
