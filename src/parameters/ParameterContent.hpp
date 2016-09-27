#ifndef _PARAMETER_CONTENT
#define _PARAMETER_CONTENT

#include "BranchPlain.hpp"
#include "BranchLength.hpp"
#include "ProtModel.hpp"
#include "Serializable.hpp"

#include <iostream>
#include <vector>

// for really properly implementing this, we'd finally need a type for
// internal and absolute branch lengths (different types of branches).
// Until then: I know, it's suboptimal

///////////////////////////////////////////////////////////////////////////////
//                             PARAMETER CONTENT                             //
///////////////////////////////////////////////////////////////////////////////
class ParameterContent : public Serializable
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    // public stuff that should be private
    std::vector<double>       values;
    std::vector<BranchPlain>  topology;
    std::vector<BranchLength> branchLengths;
    std::vector<ProtModel>    protModel;
    std::vector<nat>          rateAssignments;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ParameterContent(
        std::vector<double>      valuesI = std::vector<double>{},
        std::vector<BranchPlain> topoI = std::vector<BranchPlain>{},
        std::vector<BranchLength>blI = std::vector<BranchLength>{},
        std::vector<ProtModel>   pmI = std::vector<ProtModel>{},
        std::vector<nat>         rA = std::vector<nat>{});
    // ________________________________________________________________________
    virtual ~ParameterContent(){}
    // ________________________________________________________________________
    virtual void                           deserialize(
        std::istream&in);
    // ________________________________________________________________________
    virtual void                           serialize(
        std::ostream&out)  const;
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&          out,
        const ParameterContent&rhs);
};


#endif
