#ifndef _INSERTION_SCORE
#define _INSERTION_SCORE

///////////////////////////////////////////////////////////////////////////////
//                              INSERTION SCORE                              //
///////////////////////////////////////////////////////////////////////////////
class InsertionScore
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    InsertionScore(
        BranchPlain     _b,
        std::vector<nat>_tmp)
        : b(_b)
        , partitionParsimony(_tmp)
    {}
    // ________________________________________________________________________
    BranchPlain                            getBranch() const
    {return b;}
    // ________________________________________________________________________
    double                                 getWeight() const
    {return logProb;}
    // ________________________________________________________________________
    void                                   setWeight(
        double w)
    {logProb = w;}
    // ________________________________________________________________________
    nat                                    getScore() const
    {
        nat result = 0;

        for (auto b : partitionParsimony)
            result += b;

        return result;
    }
    // ________________________________________________________________________
    nat                                    getPartitionScore(
        int model) const {return partitionParsimony[model];}

    ///////////////////////////////////////////////////////////////////////////
    //                           PRIVATE INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ________________________________________________________________________
    friend std::ostream&                   operator<<(
        std::ostream&        out,
        const InsertionScore&rhs)
    {
        out <<  "(" << rhs.b << "=";

        for (auto elem : rhs.partitionParsimony)
            out << elem << ",";

        return out;
    }

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    BranchPlain b;
    std::vector<nat>partitionParsimony;
    double logProb;
};


#endif
