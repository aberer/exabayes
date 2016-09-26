#ifndef _BLOCK_PROPOSALCONFIG_H
#define _BLOCK_PROPOSALCONFIG_H

#include <cassert>
#include <map>

#include "ExaBlock.hpp"
#include "ProposalType.hpp"
#include "GlobalVariables.hpp"
#include "TopoMove.hpp"


// TODO allow for scientific doubles

///////////////////////////////////////////////////////////////////////////////
//                           BLOCK PROPOSAL CONFIG                           //
///////////////////////////////////////////////////////////////////////////////
class BlockProposalConfig : public ExaBlock
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BlockProposalConfig();
    // ________________________________________________________________________
    virtual void                    Read(
        NxsToken&token);
    // ________________________________________________________________________
    bool                            wasSetByUser(
        ProposalType type) const
    {return userValue.find(type) != userValue.end(); }
    // ________________________________________________________________________
    double                          getProposalWeight(
        ProposalType type) const
    {
        assert(userValue.find(type) != userValue.end());
        return userValue.at(type);
    }
    // ________________________________________________________________________
    double                          getEsprStopProp() const
    {return esprStopProp; }
    // ________________________________________________________________________
    double                          getEtbrStopProb() const
    {return etbrStopProb; }
    // ________________________________________________________________________
    double                          getParsimonyWarp() const
    {return parsimonyWarp; }
    // ________________________________________________________________________
    void                            verify();
    // ________________________________________________________________________
    int                             getParsSPRRadius() const
    {return parsSPRRadius; }
    // ________________________________________________________________________
    int                             getLikeSprMaxRadius() const
    {return _likeSprMaxRadius; }
    // ________________________________________________________________________
    double                          getLikeSprWarp() const
    {return _likeSprWarp; }
    // ________________________________________________________________________
    MoveOptMode                     getMoveOptMode()  const
    {return _moveOptMode; }
    // ________________________________________________________________________
    bool                            hasUseMultiplier() const
    {return _useMultiplier; }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::unordered_map<ProposalType, double, ProposalTypeHash> userValue;

    double etbrStopProb;
    double esprStopProp;
    double parsimonyWarp;
    int _likeSprMaxRadius;
    int parsSPRRadius;
    double _likeSprWarp;
    MoveOptMode _moveOptMode;
    bool _useMultiplier;
};


#endif
