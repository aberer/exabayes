#ifndef BRANCH_BACKUP_HPP
#define  BRANCH_BACKUP_HPP

#include "BranchPlain.hpp"
#include "BranchLength.hpp"

#include <vector>

class AbstractParameter;
class TreeAln;

///////////////////////////////////////////////////////////////////////////////
//                              BRANCH BACKUP                               //
///////////////////////////////////////////////////////////////////////////////
class BranchBackup
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    BranchBackup()
        : _backup{}
    {}
    // ________________________________________________________________________
    BranchBackup(
        TreeAln&                       traln,
        std::vector<BranchPlain>       bs,
        std::vector<AbstractParameter*>params);
    // ________________________________________________________________________
    void                                              extend(
        TreeAln&           traln,
        AbstractParameter* param,
        const BranchLength&length);
    // ________________________________________________________________________
    void                                              extend(
        TreeAln&           traln,
        AbstractParameter* param,
        const BranchPlain& length);
    // ________________________________________________________________________
    void                                              resetFromBackup(
        TreeAln& traln) const;

    // ________________________________________________________________________
    std::tuple<bool, BranchLength>                    find(
        const BranchPlain&branch,
        AbstractParameter*param) const;
    // ________________________________________________________________________
    friend std::ostream& operator                     <<(
        std::ostream&       out,
        const BranchBackup& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                             PRIVATE DATA                              //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<std::pair<AbstractParameter*, BranchLength> > _backup;
};


#endif
