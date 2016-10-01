#ifndef _BLOCK_PARTITION_H
#define _BLOCK_PARTITION_H

#include "ExaBlock.hpp"

#include "GlobalVariables.hpp"
#include "AbstractParameter.hpp"
#include "TreeAln.hpp"
#include "Category.hpp"

///////////////////////////////////////////////////////////////////////////////
//                                BLOCK PARAMS                               //
///////////////////////////////////////////////////////////////////////////////
class BlockParams : public ExaBlock
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BlockParams()
        : parameters{}
        , tralnPtr{nullptr}
    {NCL_BLOCKTYPE_ATTR_NAME = "PARAMS";}
    // ________________________________________________________________________
    BlockParams(
        const BlockParams&rhs)  = default;
    // ________________________________________________________________________
    BlockParams&                                              operator=(
        const BlockParams&rhs)  = default;
    // ________________________________________________________________________
    void                                                      setTree(
        const TreeAln* _traln){tralnPtr = _traln;}
    // ________________________________________________________________________
    vector<AbstractParameter::UPtr>                           getParameters()
    const;
    // ________________________________________________________________________
    virtual void                                              Read(
        NxsToken&token);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                                      partitionError(
        nat    partition,
        size_t totalPart) const;
    // ________________________________________________________________________
    void                                                      parseScheme(
        NxsToken& token,
        Category  cat,
        nat&      idCtr);
    // ________________________________________________________________________
    nat                                                       getNumSeen(
        Category cat);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    vector<unique_ptr<AbstractParameter> >parameters;
    const TreeAln*                         tralnPtr; // NON-owning
};


#endif
