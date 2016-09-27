#ifndef BRANCHLENGTH_H
#define BRANCHLENGTH_H

#include "BranchPlain.hpp"
#include "InternalBranchLength.hpp"
#include "Serializable.hpp"


///////////////////////////////////////////////////////////////////////////////
//                               BRANCH LENGTH                               //
///////////////////////////////////////////////////////////////////////////////
class BranchLength :  public BranchPlain
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // overwrites BranchPlain
    // ________________________________________________________________________
    //
    virtual void                            deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                            serialize(
        std::ostream&out) const;
    // ________________________________________________________________________
    explicit BranchLength(
        const BranchPlain&   b = BranchPlain(),
        InternalBranchLength len = 0)
        : BranchPlain(b)
        , _length(len)
    {}
    // ________________________________________________________________________
    virtual ~BranchLength(){}
    // ________________________________________________________________________
    double                                  toMeanSubstitutions(
        double meanSubstRate) const
    {return _length.toMeanSubstitutions(meanSubstRate); }
    // ________________________________________________________________________
    InternalBranchLength                    getLength() const
    {return _length; }
    // ________________________________________________________________________
    void                                    setLength(
        InternalBranchLength length)
    {_length = length; }
    // ________________________________________________________________________
    BranchLength                            getInverted() const
    {return BranchLength(BranchPlain::getInverted(), _length); }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    InternalBranchLength _length;
};


#endif /* BRANCHLENGTH_H */
