#ifndef _TREE_PROCESSOR_HPP
#define _TREE_PROCESSOR_HPP

#include "TreeAln.hpp"

#include <string>


///////////////////////////////////////////////////////////////////////////////
//                               TREE PROCESSOR                              //
///////////////////////////////////////////////////////////////////////////////
class TreeProcessor
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    TreeProcessor(
        std::vector<std::string>fileNames,
        bool                    expensiveCheck);
    // ________________________________________________________________________
    TreeProcessor(
        TreeProcessor&& tp);
    // ________________________________________________________________________
    virtual ~TreeProcessor(){}
    // ________________________________________________________________________
    TreeProcessor&                                    operator=(
        TreeProcessor&&tp);
    // ________________________________________________________________________
    const std::vector<std::string>                    getTaxa() const
    {return _taxa;}

    ///////////////////////////////////////////////////////////////////////////
    //                          PROTECTED INTERFACE                          //
    ///////////////////////////////////////////////////////////////////////////
protected:
    // ________________________________________________________________________
    auto                                              fillTaxaInfo(
        std::string fileName)->std::vector<std::string>;
    // ________________________________________________________________________
    template<bool readBl>
    void                                              nextTree(
        std::istream&treefile);
    // ________________________________________________________________________
    void                                              skipTree(
        std::istream&iss);
    // ________________________________________________________________________
    static std::string                                trim(
        const std::string& str,
        const std::string& whitespace  = " \t");

    ///////////////////////////////////////////////////////////////////////////
    //                             PROTECTED DATA                            //
    ///////////////////////////////////////////////////////////////////////////
protected:
    std::unique_ptr<TreeAln>_tralnPtr;
    std::vector<std::string>_fns;
    std::vector<std::string>_taxa;
};


#endif
