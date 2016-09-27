#ifndef _PARTITION_ASSIGNMENT_HPP
#define _PARTITION_ASSIGNMENT_HPP

#include "common.h"
#include "Partition.hpp"

#include <map>

///////////////////////////////////////////////////////////////////////////////
//                                 ASSIGNMENT                                //
///////////////////////////////////////////////////////////////////////////////
struct Assignment
{
    nat               partId;
    nat               procNum;
    nat               offset;
    nat               width;
    nat               compStates;
};


// ____________________________________________________________________________
std::ostream&
                                                                                                                                            operator
<<(
    std::ostream&     out,
    const Assignment& rhs);

///////////////////////////////////////////////////////////////////////////////
//                                 PART INFO                                 //
///////////////////////////////////////////////////////////////////////////////
struct PartInfo
{
    nat               id;
    nat               num;
    nat               compStates; // not used in modern thingie
};


///////////////////////////////////////////////////////////////////////////////
//                            PARTITION ASSIGNMENT                           //
///////////////////////////////////////////////////////////////////////////////
class PartitionAssignment
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    explicit PartitionAssignment(
        size_t size)
        : _numProc(size)
        , _proc2assignment{}
    {}
    // ________________________________________________________________________
    size_t
                                                                                                                                            getTotalWidthPerProc(
        nat proc) const;
    // ________________________________________________________________________
    void
                                                                                                                                            assignOld(
        const std::vector<Partition>& pass);
    // ________________________________________________________________________
    void
                                                                                                                                            assign(
        const std::vector<Partition>& pass);
    // ________________________________________________________________________
    const std::multimap<nat,
                        Assignment>&
    getAssignment()
    const
    {return _proc2assignment;}
    // ________________________________________________________________________
    auto
                                                                                                                                            getCountsAndDispls(
        size_t bla) const
        ->std::pair<std::vector<int>, std::vector<int> >;
    // ________________________________________________________________________
    size_t
                                                                                                                                            getNumProc()
    const
    {return _numProc;}
    // ________________________________________________________________________
    auto
                                                                                                                                            getNumPartPerProcess()
    const->std::vector<nat>;
    // ________________________________________________________________________
    auto
                                                                                                                                            getSitesPerProcess()
    const->std::vector<nat>;

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void
                                                                                                                                            assignForType(
        std::vector<PartInfo>);
    // ________________________________________________________________________
    void
                                                                                                                                            improvedAssign(
        std::vector<PartInfo>partitions);
    // ________________________________________________________________________
    void
                                                                                                                                            _assignToProcFull(
        PartInfo         p,
        nat              proc,
        std::vector<nat>&numAssigned,
        std::vector<nat>&sizeAssigned);
    // ________________________________________________________________________
    void
                                                                                                                                            _assignToProcPartially(
        PartInfo         p,
        nat              proc,
        std::vector<nat>&numAssigned,
        std::vector<nat>&sizeAssigned,
        nat              offset,
        nat              numElem);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    size_t _numProc;
    std::multimap<nat, Assignment>_proc2assignment;
};


#endif
