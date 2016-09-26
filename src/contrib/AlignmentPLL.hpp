#ifndef PHYLIPALIGNMENT_H
#define PHYLIPALIGNMENT_H

#include "pll.h"
#include "AbstractAlphabet.hpp"

#include <iostream>
#include <memory>

///////////////////////////////////////////////////////////////////////////////
//                               ALIGNMENT PLL                               //
///////////////////////////////////////////////////////////////////////////////
class AlignmentPLL
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    AlignmentPLL();
    // ________________________________________________________________________
    AlignmentPLL(
        const AlignmentPLL& rhs) = delete;
    // ________________________________________________________________________
    AlignmentPLL&                   operator=(
        const AlignmentPLL&rhs)  = delete;
    // ________________________________________________________________________
    AlignmentPLL(
        AlignmentPLL&& rhs);
    // ________________________________________________________________________
    AlignmentPLL&                   operator=(
        AlignmentPLL&&rhs);
    // ________________________________________________________________________
    void                            substituteBases();
    // ________________________________________________________________________
    void                            initAln(
        std::string alnFile,
        int         fileType);
    // ________________________________________________________________________
    void                            initPartitions(
        std::string partitionFile);
    // ________________________________________________________________________
    void                            print() const;
    // ________________________________________________________________________
    virtual ~AlignmentPLL();
    // ________________________________________________________________________
    void                            writeToFile(
        std::string fileName) const;
    // ________________________________________________________________________
    void                            writeHeader(
        std::ofstream&out) const;
    // ________________________________________________________________________
    template<typename T>
    void                            myWrite(
        std::ostream& out,
        const T*      ptr,
        size_t        length) const
    {out.write((const char*) ptr, sizeof(T) * length); }
    // ________________________________________________________________________
    void                            writeWeights(
        std::ofstream&out) const;
    // ________________________________________________________________________
    void                            createDummyPartition(
        Alphabet alphabet);
    // ________________________________________________________________________
    static int                      guessFormat(
        std::string file);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    pllAlignmentData* _pllAlignmentData;
    partitionList*    _partitions;
};


#endif /* PHYLIPALIGNMENT_H */
