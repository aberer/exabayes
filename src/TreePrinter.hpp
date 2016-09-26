#ifndef _TREE_PRINTER_HPP
#define _TREE_PRINTER_HPP

#include "TreeAln.hpp"

#include <sstream>

class AbstractParameter;

///////////////////////////////////////////////////////////////////////////////
//                                TREE PRINTER                               //
///////////////////////////////////////////////////////////////////////////////
class TreePrinter
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // _________________________________________________________________________
    TreePrinter(
        bool withBranchLengths,
        bool withInternalNodes,
        bool withRealNames);
    // _________________________________________________________________________
    //  prints the tree belonging to model
    //
    std::string                    printTree(
        const TreeAln&                        traln,
        const std::vector<AbstractParameter*>&params);
    // _________________________________________________________________________
    std::string                    printTree(
        const TreeAln&     traln,
        AbstractParameter* params);
    // _________________________________________________________________________
    std::string                    printTree(
        const TreeAln& traln);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // _________________________________________________________________________
    void                           helper(
        const TreeAln&                         traln,
        std::stringstream&                     ss,
        nodeptr                                p,
        bool                                   isFirst,
        const std::vector<AbstractParameter*>& params);
    // _________________________________________________________________________
    void                           printBranchLength(
        const TreeAln&                         traln,
        std::stringstream&                     ss,
        nodeptr                                p,
        const std::vector<AbstractParameter*>& params);


    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    bool withBranchLengths;
    bool withInternalNodes;
    bool withRealNames;
};


#endif
