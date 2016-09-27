#ifndef BRANCHLENGTHS_H
#define BRANCHLENGTHS_H

#include <vector>

#include "BranchPlain.hpp"
#include "InternalBranchLength.hpp"
#include "Serializable.hpp"

///////////////////////////////////////////////////////////////////////////////
//                               BRANCH LENGTHS                              //
///////////////////////////////////////////////////////////////////////////////
class BranchLengths : public BranchPlain
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit BranchLengths(
        const BranchPlain&               b = BranchPlain(),
        std::vector<InternalBranchLength>lengths = {{}})
        : BranchPlain(b)
        , _lengths(lengths)
    {}
    // ________________________________________________________________________
    virtual ~BranchLengths(){}
    // ________________________________________________________________________
    BranchLengths                                        getInverted() const
    {return BranchLengths(BranchPlain::getInverted(), _lengths); }
    // ________________________________________________________________________
    void                                                 setLengths(
        std::vector<InternalBranchLength>lengths)
    {_lengths = lengths; }
    // ________________________________________________________________________
    std::vector<InternalBranchLength>                    getLengths() const
    {return _lengths; }

    // implements Serializable
    // ________________________________________________________________________
    virtual void                                         deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                         serialize(
        std::ostream&out) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<InternalBranchLength> _lengths;
};


#endif /* BRANCHLENGTHS_H */
