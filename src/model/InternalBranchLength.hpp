#ifndef INTERNALBRANCHLENGTH_H
#define INTERNALBRANCHLENGTH_H

#include "Serializable.hpp"

///////////////////////////////////////////////////////////////////////////////
//                           INTERNAL BRANCH LENGTH                          //
///////////////////////////////////////////////////////////////////////////////
class InternalBranchLength : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    InternalBranchLength(
        double z = 0)
        : _internalBranchLength{z}
    {}

    // implements Serializable
    // ________________________________________________________________________
    virtual void                                   deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                                   serialize(
        std::ostream&out) const;

    // ________________________________________________________________________
    virtual ~InternalBranchLength(){}

    // ________________________________________________________________________
    /**
     *  Creates an internal branch length (i.e., z-value) for a tree
     */
    static InternalBranchLength                    fromAbsolute(
        double absLen,
        double meanSubstRate);

    // ________________________________________________________________________
    double                                         getValue() const
    {return _internalBranchLength; }

    // ________________________________________________________________________
    void                                           setValue(
        double val)
    {_internalBranchLength = val; }

    // ________________________________________________________________________
    double                                         toMeanSubstitutions(
        double meanSubstRate) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    double _internalBranchLength;
};


#endif /* INTERNALBRANCHLENGTH_H */
