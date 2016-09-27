#ifndef _BASIC_TREE_READER
#define _BASIC_TREE_READER

#include "BranchLength.hpp"

#include "BranchLengthPolicy.hpp"
#include "LabelPolicy.hpp"
#include <iosfwd>

#include <string>
#include <vector>
#include <sstream>


typedef unsigned int nat;

///////////////////////////////////////////////////////////////////////////////
//                             BASIC TREE READER                             //
///////////////////////////////////////////////////////////////////////////////
template<class LABEL_READER, class BL_READER>
class BasicTreeReader
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BasicTreeReader(
        nat numTax);
    // ________________________________________________________________________
    std::vector<BranchLength>                    extractBranches(
        std::istream&iss);
    // ________________________________________________________________________
    void                                         setLabelMap(
        std::unordered_map<std::string, nat>map)
    {lr.setLabelMap(map);}

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    double                                       parseFloat(
        std::istream&iss);
    // ________________________________________________________________________
    std::tuple<nat, double>                      parseSubTree(
        std::istream&             iss,
        std::vector<BranchLength>&branches,
        bool);
    // ________________________________________________________________________
    std::tuple<nat, double>                      parseElement(
        std::istream&iss);
    // ________________________________________________________________________
    void                                         expectChar(
        std::istream&iss,
        int          ch);
    // ________________________________________________________________________
    void                                         addBranch(
        nat label,
        std::tuple<nat, double>subtree,
        std::vector<BranchLength>&branches) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    nat _highestInner;
    LABEL_READER lr;
    BL_READER br;
};


#endif
