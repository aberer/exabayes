#ifndef _BYTE_FILE_HPP
#define _BYTE_FILE_HPP

#include "Partition.hpp"

#include <fstream>
#include <string>

class Communicator;
class ParallelSetup;
class PartitionAssignment;


///////////////////////////////////////////////////////////////////////////////
//                                  POSITION                                 //
///////////////////////////////////////////////////////////////////////////////
enum class Position
{
    HEADER,
    WEIGHTS,
    TAXA,
    PARTITIONS,
    ALIGNMENT,
    PARSIMONY
};


///////////////////////////////////////////////////////////////////////////////
//                               INTEGER WIDTH                               //
///////////////////////////////////////////////////////////////////////////////
enum class IntegerWidth : size_t
{
    BIT_64 = sizeof(uint64_t),
    BIT_32 = sizeof(uint32_t),
    BIT_16 = sizeof(uint16_t),
    BIT_8  = sizeof(uint8_t)
};


///////////////////////////////////////////////////////////////////////////////
//                                 BYTE FILE                                 //
///////////////////////////////////////////////////////////////////////////////
class ByteFile
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ByteFile(
        std::string name,
        bool        saveMemory = false);
    // ________________________________________________________________________
    void                                            parse(
        ParallelSetup&pl);
    // ________________________________________________________________________
    std::vector<std::string>                        getTaxa() const
    {return _taxa;}
    // ________________________________________________________________________
    auto                                            getPartitions() const
        ->std::vector<Partition>;
    // ________________________________________________________________________
    nat                                             determineOptStride(
        Communicator& comm);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    auto                                            parseHeader(
        ParallelSetup& pl)
        ->std::tuple<int, int>;
    // ________________________________________________________________________
    void                                            checkMagicNumber();
    // ________________________________________________________________________
    void                                            parsePartitions(
        int numPart);
    // ________________________________________________________________________
    void                                            seek(
        Position pos);
    // ________________________________________________________________________
    void                                            parseTaxa(
        int numTax);
    // ________________________________________________________________________
    template<typename T>
    void                                            parseWeights(
        ParallelSetup&       pl,
        PartitionAssignment& pAss);
    // ________________________________________________________________________
    void                                            parseAlns(
        ParallelSetup&       pl,
        PartitionAssignment& pAss);
    // ________________________________________________________________________
    void                                            distributeArray(
        ParallelSetup& pl);
    // ________________________________________________________________________
    void                                            parseAlnsDirectRead(
        ParallelSetup&       pl,
        PartitionAssignment& pAss);
    // ________________________________________________________________________
    void                                            parseAlnsDirect_newLayout(
        ParallelSetup&       pl,
        PartitionAssignment& pAss);
    // ________________________________________________________________________
    template<typename T>
    std::vector<std::vector<T> >                    readAndDistributeArrays(
        ParallelSetup&       pl,
        PartitionAssignment& pAss,
        size_t               numberOfArrays);
    // ________________________________________________________________________
    template<typename T>
    T                                               readVar();
    // ________________________________________________________________________
    template<typename T>
    std::vector<T>                                  readArray(
        size_t length);
    // ________________________________________________________________________
    template<typename T>
    void                                            readArray(
        size_t length,
        T*     array);

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::string _fileName;
    std::ifstream _in;
    std::vector<std::string>_taxa;
    int _numPat;
    std::vector<Partition>_partitions;
    bool _saveMemory;
    IntegerWidth _weightType;
    // CLOCK::system_clock::time_point _initTime;
};


#include "ByteFile.tpp"

#endif
