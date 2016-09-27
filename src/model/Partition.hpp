#ifndef _PARTITION_HPP
#define _PARTITION_HPP

#include "pll.h"

#include "extensions.hpp"

#include <memory>
#include <iosfwd>

///////////////////////////////////////////////////////////////////////////////
//                                 PARTITION                                 //
///////////////////////////////////////////////////////////////////////////////
class Partition
{
    ///////////////////////////////////////////////////////////////////////////
    //                              PUBLIC DATA                              //
    ///////////////////////////////////////////////////////////////////////////
public:
    static int           maxCategories; // TODO
    static constexpr int INTS_PER_VECTOR = size_t(EXA_ALIGN) / sizeof(nat);

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    Partition(
        nat         numTax,
        std::string name,
        int         dataType,
        int         states,
        int         maxTipStates,
        bool        saveMemory);
    // ________________________________________________________________________
    ~Partition(){}
    // ________________________________________________________________________
    Partition(
        const Partition&rhs);
    // ________________________________________________________________________
    Partition(
        Partition&&rhs) = default;
    // ________________________________________________________________________
    Partition&                             operator=(
        Partition rhs);
    // ________________________________________________________________________
    friend void                            swap(
        Partition& lhs,
        Partition& rhs);
    // ________________________________________________________________________
    void                                   defaultInit();
    // ________________________________________________________________________
    pInfo&                                 getHandle()
    {return _partition; }
    // ________________________________________________________________________
    const pInfo&                           getHandle() const
    {return _partition;    }
    // ________________________________________________________________________
    std::string                            getName() const
    {return _name; }
    // ________________________________________________________________________
    int                                    getStates() const
    {return _partition.states; }
    // ________________________________________________________________________
    void                                   setStates(
        int states){_partition.states = states; }
    // ________________________________________________________________________
    int                                    getMaxTipStates() const
    {return _partition.maxTipStates; }
    // ________________________________________________________________________
    void                                   setMaxTipStates(
        int s)
    {_partition.maxTipStates = s; }
    // ________________________________________________________________________
    int                                    getLower() const
    {return _partition.lower; }
    // ________________________________________________________________________
    void                                   setLower(
        int lower)
    {_partition.lower = lower; }
    // ________________________________________________________________________
    int                                    getUpper() const
    {return _partition.upper; }
    // ________________________________________________________________________
    void                                   setUpper(
        int upper)
    {_partition.upper = upper; }
    // ________________________________________________________________________
    int                                    getWidth() const
    {return _partition.width; }
    // ________________________________________________________________________
    void                                   setWidth(
        int width)
    {_partition.width = width; }
    // ________________________________________________________________________
    int                                    getDataType() const
    {return _partition.dataType; }
    // ________________________________________________________________________
    void                                   setDataType(
        int dataType)
    {_partition.dataType = dataType; }
    // ________________________________________________________________________
    int                                    getProtModels() const
    {return _partition.protModels; }
    // ________________________________________________________________________
    void                                   setProtModel(
        int protModel)
    {_partition.protModels = protModel; }
    // ________________________________________________________________________
    int                                    getProtFreqs() const
    {return _partition.protUseEmpiricalFreqs; }
    // ________________________________________________________________________
    void                                   setProtFreqs(
        int protFreqs)
    {_partition.protUseEmpiricalFreqs = protFreqs; }
    // ________________________________________________________________________
    void                                   setNonGTR(
        int nonGTR)
    {_partition.nonGTR = nonGTR; }
    // ________________________________________________________________________
    int                                    getNonGTR() const
    {return _partition.nonGTR; }
    // ________________________________________________________________________
    double                                 getPartitionContribution() const
    {return _partition.partitionContribution; }
    // ________________________________________________________________________
    void                                   setPartitionContribution(
        double tmp){_partition.partitionContribution = tmp; }
    // ________________________________________________________________________
    void                                   setAlignment(
        shared_pod_ptr<unsigned char>aln,
        int                          width);
    // ________________________________________________________________________
    void                                   setWeights(
        shared_pod_ptr<int>wgts,
        int                width);
    // ________________________________________________________________________
    double                                 getAlpha() const
    {return _partition.alpha; }
    // ________________________________________________________________________
    void                                   setAlpha(
        double a)
    {_partition.alpha = a; }
    // ________________________________________________________________________
    std::vector<double>                    getSubstRates() const
    {return _substRates; }
    // ________________________________________________________________________
    void                                   setSubstRates(
        std::vector<double>r)
    {_substRates = r; }
    // ________________________________________________________________________
    std::vector<double>                    getFrequencies() const
    {return _frequencies; }
    // ________________________________________________________________________
    void                                   setFrequencies(
        std::vector<double>f)
    {_frequencies = f; }
    // ________________________________________________________________________
    void                                   prepareParsimony();
    // ________________________________________________________________________
    double                                 getFracChange() const
    {return _partition.fracchange; }
    // ________________________________________________________________________
    void                                   setFracChange(
        double f)
    {_partition.fracchange = f; }
    // ________________________________________________________________________
    nat                                    getPartitionParsimony()
    {return _parsimonyScore.at(0); }
    // ________________________________________________________________________
    std::ostream&                          printAlignment(
        std::ostream&out);
    // ________________________________________________________________________
    int                                    getGapVectorLength() const
    {return _partition.gapVectorLength; }
    // ________________________________________________________________________
    void                                   setGapVectorLength(
        int len)
    {_partition.gapVectorLength = len; }
    // ________________________________________________________________________
    void                                   setParsimonyInformative(
        std::vector<bool>theInfo)
    {_parsimonyInformative = theInfo; }
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&    out,
        const Partition& rhs);

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    void                                   setPtrs();
    // ________________________________________________________________________
    bool                                   isInformative(
        int site) const;
    // ________________________________________________________________________
    std::vector<bool>                      determineUninformativeSites() const;
    // ________________________________________________________________________
    void                                   compressDNA(
        const std::vector<bool>&informative);
    // ________________________________________________________________________
    void                                   initGapVectorStruct();

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    bool                          _saveMemory;
    nat                           _numTax;
    pInfo                         _partition;
    std::string                   _name;

    shared_pod_ptr<int>           _wgtPtr;
    shared_pod_ptr<unsigned char> _y;

    typename aligned_vector<double>::type _left;
    typename aligned_vector<double>::type _right;
    typename aligned_vector<double>::type _EV;
    typename aligned_vector<double>::type _tipVector;
    std::vector<double>           _EIGN;
    std::vector<double>           _EI;
    std::vector<double>           _substRates;
    std::vector<double>           _frequencies;
    std::vector<double>           _empiricalFrequencies;
    std::vector<double>           _gammaRates;
    std::vector<nat>              _globalScaler;
    std::vector<unsigned char*>   _yPtrs;

    std::vector<double*>          _xVector;
    std::vector<size_t>           _xSpaceVector;

    shared_pod_ptr<nat>           _parsVect;
    std::vector<nat>              _parsimonyScore;

    std::vector<bool>             _parsimonyInformative;
    std::vector<nat>              _gapVector;
    typename aligned_vector<double>::type _gapColumn;
    typename aligned_vector<double>::type _sumBuffer;
};


#endif
